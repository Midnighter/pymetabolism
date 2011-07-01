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


import warnings

from .metabolism import metabolism


class SBMLParser(object):
    """A class implementing methods for parsing a SBML model

    @ivar path: path to SBML model
    @type path: C{str}

    @ivar rprefix: prefix e.g. 'R_' that should be stripped from reaction identifiers
    @type rprefix: C{str} or compiled regex object e.g. re.compile('R.')

    @ivar rsuffix: suffix e.g. '_Ec_Coli_Core' that should be stripped from reaction identifiers
    @type rsuffix: C{str} or compiled regex object e.g. re.compile('_.*$')

    @ivar cprefix: prefix e.g. 'M_' that should be stripped from compound identifiers
    @type cprefix: C{str} or compiled regex object e.g. re.compile('M.')

    @ivar csuffix: suffix e.g. '_c' that should be stripped from compound identifiers
    @type csuffix: C{str} or compiled regex object e.g. re.compile('_.$')

    @todo: implement convenience stuff
    @todo: Fix suffix and prefix handling
    @todo: Use logging
    @todo: Check for boundary conditions
    """

    _sbml = None

    def __init__(self, reaction_prefix="R_", reversible_suffix="r",
            compound_prefix="M_", compartment_suffix={"_c": "Cytosol", "_e":
            "Extra_organism", "_b": "Exchange", "_p": "Periplasm"}):
        if not self.__class__._sbml:
            try:
                self.__class__._sbml = __import__(name="libsbml")
            except ImportError:
                raise ImportError("libsbml is required for this functionality,"\
                        " see http://sbml.org/Software/libSBML")
        object.__init__(self)
        self.reaction_prefix = reaction_prefix
        self.reversible_suffix = reversible_suffix
        self.compound_prefix = compound_prefix
        self.compartment_suffix = compartment_suffix
        self.compartments = None
        self.compounds = None
        self.reactions = None

    def parse(self, path):
        """
        """
        document = self._sbml.readSBML(path)
        if document.getNumErrors() > 0:
            warnings.warn("reading the SBML document '%s' produced some errors"\
                    % path)
        model = document.getModel()
        # parse compartments
        self.compartments = [self._parse_compartment(comp) for comp in
                model.getListOfCompartments()]
        # parse compounds
        self.compounds = [self._parse_species(cmpd) for cmpd in
                model.getListOfSpecies()]
        self.reactions = [self._parse_reaction(rxn, model) for rxn in
                model.getListOfReactions()]

    def _parse_compartment(self, compartment):
        suffix = ""
        for (suff, name) in self.compartment_suffix.iteritems():
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
        if identifier.startswith(self.compound_prefix):
            identifier = identifier[len(self.compound_prefix):]
        compartment = None
        for suffix in self.compartment_suffix:
            if identifier.endswith(suffix):
                identifier = identifier[:-len(suffix)]
                compartment = metabolism.SBMLCompartment(self.compartment_suffix[suffix],
                        suffix=suffix)
                break
        return (identifier, compartment)

    def _parse_species(self, compound):
        """Able to parse entries from getListOfSpecies

        @todo: Check for meta information and parse if available
        """
        (identifier, comp) = self._strip_species_id(compound.getId())
        if not comp:
            comp = metabolism.SBMLCompartment(compound.getCompartment())
        cmpd = metabolism.SBMLCompound(identifier, charge=compound.getCharge())
        if not comp:
            return cmpd
        else:
            return metabolism.SBMLCompartmentCompound(cmpd, comp)

    def _strip_reaction_id(self, name):
        identifier = name
        if identifier.startswith(self.reaction_prefix):
            identifier = identifier[len(self.reaction_prefix):]
        if identifier.endswith(self.reversible_suffix):
            identifier = identifier[:-len(self.reversible_suffix)]
        return identifier

    def _parse_reaction(self, reaction, model):
        """Able to parse entries from getListOfReactions"""
        identifier = self._strip_reaction_id(reaction.getId())
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
                reversible=reaction.getReversible(), **params)

