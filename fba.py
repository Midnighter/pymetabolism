#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===========================
Flux Balance Analysis Model
===========================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
    Nikolaus Sonnenschein
:Date:
    2011-03-28
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    fba.py
"""

import numpy

class FBAModel(object):
    """
    """

    def __init__(self, model):
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
        self._model = model

    def delete_reaction(self, reaction):
        """
        Takes a list of reactions and constrains their flux to zero.

        Parameters
        ----------
        reaction: iterable
            iterable that contains the names of the reactions
        """
        reactions = dict()
        for reaction in reaction:
            reactions[reaction] = (0, 0)
        self._model.modify_column_constraints(reactions)

    def delete_reaction_from_stoich(self, reaction):
        """
        Takes a list of metabolites and removes them from the stoichiometry matrix.

        Parameters
        ----------
        reaction: iterable
            iterable that contains the names of the reactions
        """
        self._model.delete_columns(reaction)

    def delete_metabolite(self, metabolite):
        """
        Takes a list of metabolites and constrains them to zero.
        Parameters
        ----------
        metabolite: iterable
            iterable that contains the names of the compounds
        """
        metabolites = dict()
        for metabolite in metabolite:
            metabolites[metabolite] = (0, 0)
        self._model.modify_row_constraints(metabolites)

    def delete_metabolite_from_stoich(self, metabolite):
        """
        Takes a list of metabolites and removes them from the stoichiometry matrix.
        Parameters
        ----------
        name: metabolite iterable
            description: list that contains the names of the metabolites
        """
        self._model.delete_rows(metabolite)

    def free_metabolites(self, metabolite):
        """
        Takes a list of metabolites and removes their constraints.
        Parameters
        ----------
        name: metabolite iterable
            description: list that contains the names of the metabolites
        """
        boundaries = dict()
        for metabolite in metabolite:
            boundaries[metabolite] = (-numpy.inf, numpy.inf)
        self._model.modify_row_constraints(boundaries)

    def modify_reaction_bounds(self, bounds_dict):
        """
        Takes a dictionary of  boundary conditions for reaction fluxes and enforces
        the respective boundary conditions.
        Parameters
        ----------
        bounds : dictionary
            description: key = name of reaction, values = tuples of lower and
            upper bounds for the flux
        """
        self._model.modify_column_constraints(self, bounds_dict)

    def modify_metabolite_bounds(self, bounds_dict):
        """
        Takes a dictionary of  boundary conditions for metabolites and enforces
        the respective boundary conditions.
        Parameters
        ----------
        bounds : dictionary
            description: key = name of metabolite, values = tuples of lower and
            upper bounds for the flux
        """
        self._model.modify_row_constraints(self, bounds_dict)

    def add_metabolite_drains(self, metabolites, suffix="_Drain"):
        """
        Takes a list of metabolites and adds reactions that consume those metabolites.
        Parameters
        ----------
        metabolites : iterable
            description: list that contains the names of metabolites that
            need to be exported
        """
        new_columns = dict()
        coeff = dict()
        for met in metabolites:
            coeff[met][0]=-1
        for met in metabolites:
            new_columns[met+suffix] = coeff[met]
        self._model.add_columns(new_columns)

    def add_transporters(self, metabolites):
        """
        Takes a list of metabolites and adds reactions that produce those metabolites.
        Parameters
        ----------
        metabolites : iterable
            description: list that contains the names of metabolites that
            need to be transported into the cell
        """
        suffix = "_Transp"
        new_columns = dict()
        coeff = dict()
        for met in metabolites:
            coeff[met][0]=1
        for met in metabolites:
            new_columns[met+suffix] = coeff[met]
        self._model.add_columns(new_columns)

    def get_reactions(self, drain="_Drain", transp="_Transp"):
        """
        Gets the list of all reactions (excluding transport and drain reactions).

        Returns
        -------
        iterable:
            The returned list contains the names of all the reactions.
        """
        temp_react = set(self.get_column_names())
        temp_drain_transp = set()
        for react in temp_react:
            if drain in react or transp in react:
                temp_drain_transp.add(react)
        return list(temp_react.difference(temp_drain_transp))

    def get_metabolites(self):
        """
        Gets the list of all metabolites.

        Returns
        -------
        iterable:
            The returned list contains the names of all the metabolites.
        """
        return list(self.get_row_names())

    def get_transporters(self, transp="_Transp"):
        """
        Gets the list of all transporters.

        Returns
        -------
        iterable:
            The returned list contains the names of all the transport reactions.
        """
        temp_transp = set()
        for react in set(self.get_column_names()):
            if transp in react:
                temp_transp.add(react)
        return list(temp_transp)

    def get_substrates_and_products(self, reaction):
        """
        Gets the substrates and products of a reaction.

        Returns
        -------
        tuple:
            The returned tuple contains the names of the substrates and the
            names of the products of the reaction.
        """
    





















        
