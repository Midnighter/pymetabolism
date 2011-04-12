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
import pymetabolism.lpmodel as pylp

from pymetabolism.parsers import SBMLParser


class MetabolicModelBuilder(object):
    """
    """

    def __init__(self, parser="SBML", lp_solver="gurobi"):
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

    def generate_model(self, path):
        from pymetabolism.fba import FBAModel
        self._parser.parse(path)
        known_fluxes = dict()
        objectives = dict()
        for rxn in self._parser.reactions:
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
                self._model.add_column(rxn.name + "_Rev", constraints, (0.0, ub))
            if rxn.flux_value is not None:
                if rxn.flux_value > 0.0:
                    known_fluxes[rxn.name] = rxn.flux_value
                elif rxn.flux_value < 0.0:
                    known_fluxes[rxn.name + "_Rev"] = abs(rxn.flux_value)
                else:
                    known_fluxes[rxn.name] = 0.0
                    if rxn.reversible:
                        known_fluxes[rxn.name + "_Rev"] = 0.0
            if rxn.objective_coefficient:
                objectives[rxn.name] = rxn.objective_coefficient
        self._model.set_objective(objectives)
        return (FBAModel(self._model), known_fluxes)

    def generate_network(self, path):
        self._parser.parse(path)
        pass

