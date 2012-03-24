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


__all__ = ["SBMLParser"]


import os
import logging

from .metabolism import metabolism
from . import miscellaneous as misc
from .errors import PyMetabolismError


logger = logging.getLogger(__name__)
logger.addHandler(misc.NullHandler())


class SBMLParser(object):
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

    def parse(self, path):
        """
        Parse a document in SBML format.
        """
        if not os.path.exists(path):
            raise PyMetabolismError("no such file '%s'" % path)
        document = self._sbml.readSBML(path)
        if document.getNumErrors() > 0:
            logger.warn("reading the SBML document '%s' produced some errors",
                    path)
        model = document.getModel()
        # parse compartments
        compartments = [self._parse_compartment(comp) for comp in
                model.getListOfCompartments()]
        # parse compounds
        compounds = [self._parse_species(cmpd) for cmpd in
                model.getListOfSpecies()]
        reactions = [self._parse_reaction(rxn, model) for rxn in
                model.getListOfReactions()]
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

