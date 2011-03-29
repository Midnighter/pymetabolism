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

    def knock_out_reaction(self, reaction):
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

    def knock_out_metabolite(self, metabolite):
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
        metabolite: iterable
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

    def add_transporters(self, metabolites, suffix = "_Transp"):
        """
        Takes a list of metabolites and adds reactions that produce those metabolites.
        Parameters
        ----------
        metabolites : iterable
            description: list that contains the names of metabolites that
            need to be transported into the cell
        """
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
        for react in set(self._model.get_column_names()):
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
        substrates=list()
        products=list()
        coeff=self._model.get_column_coefficients(reaction)#get_column_coefficients
        # would return a dictionary with key=name of metabolite, value=coefficient
        for i, elem in coeff.items():
            if elem > 0.:
                products.append(i)
            elif elem < 0.:
                substrates.append(i)
        return (substrates,products) 
        
    def get_reaction_objective(self):
        """
        Gets the reaction objective (maximization).

        Returns
        -------
        string:
            The returned string is the reaction that is currently set for
            maximization.
        """
        objective_dict = self._model.get_objective()#the get_objective function 
        #would return a dictionary where key=name of reaction, value=0 if 
        #reaction is not objective or 1 if reaction is objective
        tmp = list()
        for k, v in objective_dict.items():
            if v == 1.:
                tmp.append(k)
        if len(tmp) is not 1:
            raise Exception, "There exists no unique reaction objective. " + str(tmp)
        return tmp[0]
        
    def set_reaction_objective(self, reaction):
        """
        Sets a certain reaction as objective (for maximization).

        Parameters:
        -------
        reaction: iterable
            iterable that contains the name of the reaction.
        """
        self._model.set_objective(reaction)
        
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
        
      
        





















