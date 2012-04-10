#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=============
Model Parsers
=============

:Authors:
    Moritz Emanuel Beber
    Nikolaus Sonnenschein
:Date:
    2011-04-07
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    parsers.py
"""


__all__ = ["parse"]


import os
import codecs
import logging

from contextlib import contextmanager
from .metabolism import metabolism
from . import miscellaneous as misc
from .errors import PyMetabolismError
from .singletonmixin import Singleton


logger = logging.getLogger(__name__)
logger.addHandler(misc.NullHandler())


class SBMLParser(Singleton):
    """
    A class implementing methods for parsing a SBML model
    """

    _sbml = None

    def __init__(self):
        if not self.__class__._sbml:
            try:
                self.__class__._sbml = __import__(name="libsbml")
            except ImportError:
                raise ImportError("libsbml is required for this functionality,"\
                        " see http://sbml.org/Software/libSBML")
        object.__init__(self)
        self._options = misc.OptionsManager.get_instance()

    def parse_string(self, xml):
        """
        Parse a document in SBML format.
        """
        document = self._sbml.readSBMLFromString(xml)
        if document.getNumErrors() > 0:
            logger.warn("reading the SBML document produced some errors")
        model = document.getModel()
        # parse compartments
        compartments = [self._parse_compartment(comp) for comp in
                model.getListOfCompartments()]
        logger.debug("approx. %d compartments", len(compartments))
        # parse compounds
        compounds = [self._parse_species(cmpd) for cmpd in
                model.getListOfSpecies()]
        logger.debug("%d compounds", len(compounds))
        reactions = [self._parse_reaction(rxn, model) for rxn in
                model.getListOfReactions()]
        logger.debug("%d reactions", len(reactions))
        return metabolism.MetabolicSystem(compartments=compartments,
                reactions=reactions, compounds=compounds)

    def _parse_compartment(self, compartment):
        suffix = ""
        for (suff, name) in self._options.compartments.iteritems():
            if name == compartment.getId():
                suffix = suff
        return metabolism.SBMLCompartment(name=compartment.getId(),
                outside=compartment.getOutside(),
                constant=compartment.getConstant(),
                spatial_dimensions=compartment.getSpatialDimensions(),
                size=compartment.getSize(), units=compartment.getUnits(),
                suffix=suffix)

    def _strip_species_id(self, name):
        identifier = name
        if identifier.startswith(self._options.compound_prefix):
            identifier = identifier[len(self._options.compound_prefix):]
        compartment = None
        for suffix in self._options.compartments.iterkeys():
            if identifier.endswith(suffix):
                identifier = identifier[:-len(suffix)]
                compartment = metabolism.SBMLCompartment(
                        self._options.compartments[suffix], suffix=suffix)
                break
        return (identifier, compartment)

    def _parse_species(self, compound):
        """
        Able to parse entries from getListOfSpecies

        @todo: Check for meta information and parse if available
        """
        (identifier, comp) = self._strip_species_id(compound.getId())
        if not comp:
            comp = metabolism.SBMLCompartment(compound.getCompartment())
        name = compound.getName()
        cmpd = metabolism.SBMLCompound(identifier, extended_name=name,
                charge=compound.getCharge())
        if not comp:
            return cmpd
        else:
            return metabolism.SBMLCompartmentCompound(cmpd, comp)

    def _strip_reaction_id(self, name):
        identifier = name
        if identifier.startswith(self._options.reaction_prefix):
            identifier = identifier[len(self._options.reaction_prefix):]
        if identifier.endswith(self._options.reversible_suffix):
            identifier = identifier[:-len(self._options.reversible_suffix)]
        return identifier

    def _parse_reaction(self, reaction, model):
        """Able to parse entries from getListOfReactions"""
        identifier = self._strip_reaction_id(reaction.getId())
        name = reaction.getName()
        params = dict()
        for param in reaction.getKineticLaw().getListOfParameters():
            params[param.getId().lower()] = param.getValue()
        substrates = dict((self._parse_species(model.getSpecies(
            elem.getSpecies())),
            abs(elem.getStoichiometry())) for elem in
            reaction.getListOfReactants())
        products = dict((self._parse_species(model.getSpecies(
            elem.getSpecies())),
            abs(elem.getStoichiometry())) for elem in
            reaction.getListOfProducts())
        return metabolism.SBMLReaction(identifier, substrates, products,
                reversible=reaction.getReversible(), extended_name=name,
                **params)


def _open_tar(path, **kw_args):
    import tarfile
    kw_args["mode"] = kw_args["mode"].strip("b")
    if isinstance(path, basestring):
        return tarfile.TarFile(name=path, mode=kw_args["mode"],
                encoding=kw_args["encoding"])
    else:
        return tarfile.TarFile(fileobj=path, mode=kw_args["mode"],
                encoding=kw_args["encoding"])

def _open_gz(path, **kw_args):
    import gzip
    if isinstance(path, basestring):
        return gzip.GzipFile(filename=path, mode=kw_args["mode"])
    else:
        return gzip.GzipFile(fileobj=path, mode=kw_args["mode"])

def _open_bz2(path, **kw_args):
    import bz2
    return bz2.BZ2File(path)

def _open_zip(path, **kw_args):
    import zipfile
    kw_args["mode"] = kw_args["mode"].strip("b")
    return zipfile.ZipFile(path, mode=kw_args["mode"])

def _open_file(path, **kw_args):
    if isinstance(path, basestring):
        return codecs.open(path, mode=kw_args["mode"],
                encoding=kw_args["encoding"])
    else:
        reader = codecs.getreader(kw_args["encoding"])
        return reader(path)

archives = {"gz": _open_gz,
        "gzip": _open_gz,
        "bz2": _open_bz2
#        "zip": _open_zip,
#        "tar": _open_tar
        }


@contextmanager
def open_file(filename, **kw_args):
    path = filename
    filename = os.path.basename(filename)
    extns = filename.split(".")
    del extns[0]
    extns.reverse()
    for ext in extns:
        ext = ext.lower()
        func = archives.get(ext, _open_file)
        path = func(path, **kw_args)
    yield (path, ext)
    if not path.closed:
        path.close()


parsers = {"xml": SBMLParser,
        "sbml": SBMLParser
        }


def parse(filename, frmt=False, mode="rb", encoding="utf-8", **kw_args):
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with  open_file(filename, **kw_args) as (file_h, ext):
        if frmt:
            ext = frmt.lower()
        if ext in parsers:
            parser = parsers[ext].get_instance()
        else:
            raise PyMetabolismError("unknown metabolic system format '{}'", ext)
        system = parser.parse_string(str(file_h.read(-1)))
    return system

