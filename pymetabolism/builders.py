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


import numpy
import warnings
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
        model = self._model.copy()

#        def parse_exchange(reaction):
#            reaction.name = reaction.name[3:-3]
#            constraints = dict()
#            comp = pymet.SBMLCompartment("Exchange")
#            for cmpd in reaction.substrates:
#                if cmpd.compartment == comp:
#                    continue
#                constraints[cmpd.name] = abs(reaction.stoichiometric_coeff(cmpd))
#            for cmpd in reaction.products:
#                if cmpd.compartment == comp:
#                    continue
#                constraints[cmpd.name] = reaction.stoichiometric_coeff(cmpd)
#            if rxn.lower_bound is not None:
#                lb = rxn.lower_bound
#            else:
#                lb = -numpy.inf
#            if reaction.upper_bound is not None:
#                ub = reaction.upper_bound
#            else:
#                ub = numpy.inf
#            # exchange reactions are inversed
#            if lb >= 0.0:
#                model.add_column(rxn.name + "_Drain", constraints, (lb, ub))
#            else:
#                model.add_column(rxn.name + "_Drain", constraints, (0.0, ub))
#            for (cmpd, factor) in constraints.iteritems():
#                constraints[cmpd] = -factor
#            if lb < 0.0:
#                model.add_column(rxn.name + "_Transp", constraints, (0.0, abs(lb)))
#            else:
#                model.add_column(rxn.name + "_Transp", constraints, (0.0, ub))
#            if reaction.flux_value is not None:
#                if reaction.flux_value > 0.0:
#                    known_fluxes[reaction.name + "_Drain"] = reaction.flux_value
#                elif reaction.flux_value < 0.0:
#                    known_fluxes[reaction.name + "_Transp"] = abs(reaction.flux_value)
#                else:
#                    known_fluxes[reaction.name + "_Transp"] = 0.0
#                    known_fluxes[reaction.name + "_Drain"] = 0.0

        self._parser.parse(path)
        known_fluxes = dict()
        objectives = dict()
        for rxn in self._parser.reactions:
#            if rxn.name.startswith("EX"):
#                parse_exchange(rxn)
#                continue
            constraints = dict()
            for cmpd in rxn.substrates:
                constraints[cmpd.name] = rxn.stoichiometric_coeff(cmpd)
            for cmpd in rxn.products:
                constraints[cmpd.name] = rxn.stoichiometric_coeff(cmpd)
            if rxn.lower_bound is not None:
                lb = rxn.lower_bound
            else:
                lb = -numpy.inf
            if rxn.upper_bound is not None:
                ub = rxn.upper_bound
            else:
                ub = numpy.inf
            if lb >= 0.0:
                model.add_column(rxn.name, constraints, (lb, ub))
            else:
                model.add_column(rxn.name, constraints, (0.0, ub))
            if rxn.reversible:
                for (cmpd, factor) in constraints.iteritems():
                    constraints[cmpd] = -factor
                if abs(lb) == 0.0:
#                    warnings.warn("LB is 0.0 for reversible reaction %s, setting to UB (=%e)" % (rxn.name, ub))
                    rev_ub = ub
                else:
                    rev_ub = abs(lb)
                    
                model.add_column(rxn.name + self.suffix, constraints, (0.0,
                        rev_ub))
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
        model.set_objective(objectives)
        return (FBAModel(model), known_fluxes)

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

