#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=======================
Metabolic Model Builder
=======================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
    Nikolaus Sonnenschein
:Date:
    2011-04-03
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    builder.py
"""


import copy
import numpy
import pymetabolism.lpmodels as pylp
import pymetabolism.metabolism as pymet

from pymetabolism.parsers import SBMLParser


class MetabolicModelBuilder(object):
    """
    """

    def __init__(self, parser="SBML", lp_solver="gurobi",
            reversible_suffix="_Rev"):
        """
        description

        Parameters
        ----------
        name: type:
            description

        Returns
        -------
        type:
            description

        Notes
        -----
        text
        """
        if lp_solver.lower() == "gurobi":
            self._model = pylp.GurobiLPModelFacade()
        else:
            pass
        if parser.lower() == "sbml":
            self._parser = SBMLParser()
        self.suffix = reversible_suffix

    def generate_model(self, path):
        from pymetabolism.fba import FBAModel

        def parse_exchange(reaction):
            reaction.name = reaction.name[3:-3]
            constraints = dict()
            comp = pymet.SBMLCompartment("Exchange")
            for cmpd in reaction.substrates:
                if cmpd.compartment == comp:
                    continue
                constraints[cmpd.name] = abs(reaction.stoichiometric_coeff(cmpd))
            for cmpd in reaction.products:
                if cmpd.compartment == comp:
                    continue
                constraints[cmpd.name] = reaction.stoichiometric_coeff(cmpd)
            if reaction.upper_bound is not None:
                ub = reaction.upper_bound
            else:
                ub = numpy.inf
            self._model.add_column(reaction.name + "_Transp", constraints, (0.0, ub))
            for (cmpd, factor) in constraints.iteritems():
                constraints[cmpd] = -factor
            self._model.add_column(reaction.name + "_Drain", constraints, (0.0, ub))
            if reaction.flux_value is not None:
                if reaction.flux_value > 0.0:
                    known_fluxes[reaction.name + "_Transp"] = reaction.flux_value
                elif reaction.flux_value < 0.0:
                    known_fluxes[reaction.name + "_Drain"] = abs(reaction.flux_value)
                else:
                    known_fluxes[reaction.name + "_Transp"] = 0.0
                    known_fluxes[reaction.name + "_Drain"] = 0.0

        self._parser.parse(path)
        known_fluxes = dict()
        objectives = dict()
        for rxn in self._parser.reactions:
            if rxn.name.startswith("EX"):
                parse_exchange(rxn)
                continue
            constraints = dict()
            for cmpd in rxn.substrates:
                constraints[cmpd.name] = rxn.stoichiometric_coeff(cmpd)
            for cmpd in rxn.products:
                constraints[cmpd.name] = rxn.stoichiometric_coeff(cmpd)
            if rxn.upper_bound is not None:
                ub = rxn.upper_bound
            else:
                ub = numpy.inf
            self._model.add_column(rxn.name, constraints, (0.0, ub))
            if rxn.reversible:
                for (cmpd, factor) in constraints.iteritems():
                    constraints[cmpd] = -factor
                self._model.add_column(rxn.name + self.suffix, constraints, (0.0, ub))
            if rxn.flux_value is not None:
                if rxn.flux_value > 0.0:
                    known_fluxes[rxn.name] = rxn.flux_value
                elif rxn.flux_value < 0.0:
                    known_fluxes[rxn.name + self.suffix] = abs(rxn.flux_value)
                else:
                    known_fluxes[rxn.name] = 0.0
                    if rxn.reversible:
                        known_fluxes[rxn.name + self.suffix] = 0.0
            if rxn.objective_coefficient:
                objectives[rxn.name] = rxn.objective_coefficient
        self._model.set_objective(objectives)
        return (FBAModel(self._model), known_fluxes)

    def generate_network(self, path, name="", disjoint_reversible=False,
            stoichiometric_factors=False):
        from pymetabolism.networks import MetabolicNetwork
        self._parser.parse(path)
        net = MetabolicNetwork(name)
        for rxn in self._parser.reactions:
            for cmpd in rxn.substrates:
                if stoichiometric_factors:
                    net.add_edge(pymet.BasicCompound(cmpd.name),
                            pymet.BasicReaction(rxn.name, rxn.reversible),
                            stoichiometry=abs(rxn.stoichiometric_coeff(cmpd)))
                else:
                    net.add_edge(pymet.BasicCompound(cmpd.name),
                            pymet.BasicReaction(rxn.name, rxn.reversible))
            for cmpd in rxn.products:
                if stoichiometric_factors:
                    net.add_edge(pymet.BasicReaction(rxn.name, rxn.reversible),
                            pymet.BasicCompound(cmpd.name),
                            stoichiometry=abs(rxn.stoichiometric_coeff(cmpd)))
                else:
                    net.add_edge(pymet.BasicReaction(rxn.name, rxn.reversible),
                            pymet.BasicCompound(cmpd.name))
            if disjoint_reversible and rxn.reversible:
                for cmpd in rxn.substrates:
                    if stoichiometric_factors:
                        net.add_edge(pymet.BasicReaction(rxn.name + self.suffix,
                                rxn.reversible), pymet.BasicCompound(cmpd.name),
                                stoichiometry=abs(rxn.stoichiometric_coeff(cmpd)))
                    else:
                        net.add_edge(pymet.BasicReaction(rxn.name + self.suffix,
                                rxn.reversible), pymet.BasicCompound(cmpd.name))
                for cmpd in rxn.products:
                    if stoichiometric_factors:
                        net.add_edge(pymet.BasicCompound(cmpd.name),
                                pymet.BasicReaction(rxn.name + self.suffix,
                                rxn.reversible),
                                stoichiometry=abs(rxn.stoichiometric_coeff(cmpd)))
                    else:
                        net.add_edge(pymet.BasicCompound(cmpd.name),
                                pymet.BasicReaction(rxn.name, rxn.reversible))
        return net

