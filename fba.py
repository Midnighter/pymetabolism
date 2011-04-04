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

    def knockout_reaction(self, reaction):
        """
        Constrains the allowed reaction flux to zero.

        Parameters
        ----------
        reaction: iterable or str
            A single name or iterable that contains the names of the reactions.
        """
        if hasattr(reaction, "__iter__"):
            bounds = dict((rxn, (0.0, 0.0)) for rxn in reaction)
        else:
            bounds = {reaction: (0.0, 0.0)}
        self._model.modify_column_bounds(bounds)

    def delete_reaction_from_stoich(self, reaction):
        """
        Completely removes the reaction(s) from the model.

        Parameters
        ----------
        reaction: iterable or str
            A single name or iterable that contains the names of the reactions.
        """
        self._model.delete_column(reaction)

    def knockout_compound(self, compound):
        """
        Constrains the reaction fluxes consuming or producing elements in
        compound to zero.

        Parameters
        ----------
        compound: iterable or str
            A single name or iterable that contains the names of the compounds.
        """
        if hasattr(compound, "__iter__"):
            bounds = dict((name, (0.0, 0.0)) for cmpd in compound for name in
                    self._model.get_column_names(compound))
        else:
            bounds = dict((name, (0.0, 0.0)) for name in
                    self._model.get_column_names(compound))
        self._model.modify_column_bounds(metabolites)

    def delete_compound_from_stoich(self, compound):
        """
        Completely removes the compound(s) from the model.

        Parameters
        ----------
        reaction: iterable or str
            A single name or iterable that contains the names of the compounds.
        """
        self._model.delete_row(compound)

    def free_metabolites(self, compound):
        """
        Completely removes any bounds on compound(s).

        Parameters
        ----------
        compound: iterable or str
            A single name or iterable that contains the names of the compounds.
        """
        if hasattr(compound, "__iter__"):
            bounds = dict((name, (-numpy.inf, numpy.inf)) for name in compound)
        else:
            bounds = {compound: (-numpy.inf, numpy.inf)}
        self._model.modify_row_bounds(boundaries)

    def modify_reaction_bounds(self, bounds):
        """
        Modifies the allowed flux on reactions.

        Parameters
        ----------
        bounds : dict
            Map of reaction names to pairs with their lower and upper bound.
        """
        self._model.modify_column_bounds(self, bounds)

    def modify_compound_bounds(self, bounds):
        """
        Modifies the bounds on compounds.

        Parameters
        ----------
        bounds : dict
            Map of compound names to pairs with their lower and upper bound.
        """
        self._model.modify_row_bounds(self, bounds)

    def add_compound_drain(self, compound, suffix="_Drain", bounds=tuple()):
        """
        Adds an export channel for each compound.

        Parameters
        ----------
        compound: iterable or str
            A single name or iterable that contains the names of the compounds
            for which export from the model is added.
        suffix: str (optional)
            Identifier added to the name of drains to distinguish them from
            normal reactions.
        bounds: tuple (optional)
            A pair of lower and upper bound on the flux for all drains added.
        """
        columns = ((name + suffix, {name: -1}, bounds) for name in compound)
        self._model.add_columns(columns)

    def add_compound_transport(self, compound, suffix = "_Transp", bounds=tuple()):
        """
        Adds an export channel for each compound.

        Parameters
        ----------
        compound : iterable or str
            A single name or iterable that contains the names of the compounds
            for which export from the model is added.
        suffix: str (optional)
            Identifier added to the name of transporters to distinguish them from
            normal reactions.
        bounds: tuple (optional)
            A pair of lower and upper bound on the flux for all transporters added.
        """
        columns = ((name + suffix, {name: 1}, bounds) for name in compound)
        self._model.add_columns(columns)

    def get_reactions(self, drain="_Drain", transp="_Transp"):
        """
        Returns
        -------
        iterator:
            Iterator over all names of reactions (excluding transport and drain
            reactions).
        """
        return (rxn for rxn in self._model.get_column_names() if not
                (rxn.endswith(drain) or rxn.endswith(transp)))

    def get_compounds(self):
        """
        Returns
        -------
        iterator:
            Iterator over all names of compounds.
        """
        return self._model.get_row_names()

    def get_drains(self, drain="_Drain"):
        """
        Gets the list of all transporters.

        Returns
        -------
        iterator:
            Iterator over all names of drains.
        """
        return (rxn for rxn in self._model.get_column_names() if
                rxn.endswith(drain))

    def get_transporters(self, transp="_Transp"):
        """
        Gets the list of all transporters.

        Returns
        -------
        iterator:
            Iterator over all names of transporters.
        """
        return (rxn for rxn in self._model.get_column_names() if
                rxn.endswith(transp))

    def get_substrates_and_products(self, reaction):
        """
        Returns
        -------
        tuple:
            Pair of the names of the substrates and the names of the products
            of the reaction as lists.
        """
        substrates=list()
        products=list()
        for (cmpd, factor) in self._model.get_row_names(reaction, True):
            if factor > 0.0:
                products.append(cmpd)
            elif factor < 0.0:
                substrates.append(cmpd)
        return (substrates,products)

    def get_objective_reaction(self, coefficient=False):
        """
        Parameters
        ----------
        coefficient: bool (optional)
            Causes the returned iterator to run over pairs of reaction name and
            absolute weight in the objective function.

        Returns
        -------
        iterator:
            Current reaction(s) that are used as objectives in LP.
        """
        return self._model.get_objective(coefficient)

    def set_objective_reaction(self, reaction):
        """
        Sets a certain reaction as objective (for maximization).

        Parameters:
        -------
        reaction: dict
            Map from reaction name(s) to their factor(s).
        """
        self._model.set_objective(reaction)

    def get_objective_value(self):
        """
        Returns
        -------
        float:
            Flux of the set objective reaction(s).
        """
        growth = self._model.get_objective_value()
        return growth if growth else 0.0

    def set_reaction_objective_minimize_rest(self, reaction, factor):
        """
        Sets a certain reaction as objective (for maximization) and minimizes
        the usage of all the other reactions. All the transporter reactions are
        set to zero flux.

        Parameters:
        -------
        reaction: iterable
            description: iterable that contains the name of the reaction.
        factor: iterable
            description: iterable that contains the factor to which the fluxes 
            of the reactions that are not objective should be minimized
        """
        react_dict=dict()
        react_list=self._model.get_column_names()
        transp_list=self.get_transporters()
        for r in react_list:
            if r in transp_list:
                react_dict[r]=(0,0)
            else: 
                react_dict[r]=(factor,factor)
        self.modify_reaction_bounds(react_dict)
        self.set_reaction_objective(reaction)

    def fba(self)
        
      
        





















