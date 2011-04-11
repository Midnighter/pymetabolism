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
        columns = list()
        for rxn in self._parser.reactions:
            constraints = dict()
            for cmpd in rxn.substrates:
                constraints[cmpd.name] = rxn.stoichiometric_coeff(cmpd)
            for cmpd in rxn.products:
                constraints[cmpd.name] = rxn.stoichiometric_coeff(cmpd)
            columns.append((rxn.name, constraints, tuple()))
        self._model.add_columns(columns)
        return FBAModel(self._model)

    def generate_network(self, path):
        self._parser.parse(path)
        pass

