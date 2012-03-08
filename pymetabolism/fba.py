#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===========================
Flux Balance Analysis Model
===========================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
    Nils Kölling
    Nikolaus Sonnenschein
:Date:
    2011-03-28
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    fba.py
"""


import random
import numpy

from .lpmodels import MetaLPModelFacade


iAF1260_minimal_medium = ["ca2_b_Transp",
        "cl_b_Transp",
        "co2_b_Transp",
        "cobalt2_b_Transp",
        "cu2_b_Transp",
        "fe2_b_Transp",
        "fe3_b_Transp",
        "h_b_Transp",
        "h2o_b_Transp",
        "k_b_Transp",
        "mg2_b_Transp",
        "mn2_b_Transp",
        "mobd_b_Transp",
        "na1_b_Transp",
        "nh4_b_Transp",
        "pi_b_Transp",
        "so4_b_Transp",
        "tungs_b_Transp",
        "zn2_b_Transp",
        "cbl1_b_Transp"]


class FBAModel(object):
    """
    """

    __metaclass__ = MetaLPModelFacade
    _error = NotImplementedError("method not available for the chosen solver")

    def __copy__(self):
        raise self._error

    def __deepcopy__(self, memo=dict()):
        raise self._error

    def copy(self):
        raise self._error

    def add_reaction(self, reaction, coefficients, lb=None, ub=None):
        """
        Add a new reaction to the model.

        The method can also be used to add multiple reactions at the same time.
        In that case each argument must be an iterable of equal length and
        matching positions corresponding to the same reaction.

        Parameters
        ----------
        reaction: `BasicReaction`
            An instance of `BasicReaction`.
        coefficients: iterable
            Pairs of `BasicCompound` instances and their
            stoichiometric coefficient.
        lb: float (optional)
            A lower bound on the mass flux through this reaction. A default
            value can be set in the options.
        ub: float (optional)
            An upper bound on the mass flux through this reaction. A default
            value can be set in the options.
        """
        raise self._error

    def knockout_reaction(self, reaction):
        """
        Constrains the allowed reaction flux to zero mimicking a dysfunctional
        reaction or reactions.

        Parameters
        ----------
        reaction: iterable or `BasicReaction`
            A single `BasicReaction` instance or an iterable with multiple ones.
        """
        raise self._error

    def free_reaction(self, reaction):
        """
        Completely removes any bounds on the reaction(s).

        Parameters
        ----------
        reaction: iterable or `BasicReaction`
            A single `BasicReaction` instance or an iterable with multiple ones.
        """
        raise self._error

    def delete_reaction(self, reaction):
        """
        Completely removes the reaction(s) from the model.

        Parameters
        ----------
        reaction: iterable or `BasicReaction`
            A single `BasicReaction` instance or an iterable with multiple ones.
        """
        raise self._error

    def knockout_compound(self, compound):
        """
        Knocks out the reactions consuming or producing the given compound(s).

        Parameters
        ----------
        compound: iterable or `BasicCompound`
            A single `BasicCompound` instance or an iterable with multiple ones.
        """
        raise self._error

    def delete_compound(self, compound):
        """
        Completely removes the compound(s) from the model.

        Parameters
        ----------
        compound: iterable or `BasicCompound`
            A single `BasicCompound` instance or an iterable with multiple ones.
        """
        raise self._error

    def free_compound(self, compound):
        """
        Completely removes any bounds on the compound(s) effectively removing it
        from the model.

        Parameters
        ----------
        compound: iterable or `BasicCompound`
            A single `BasicCompound` instance or an iterable with multiple ones.
        """
        raise self._error

    def modify_reaction_bounds(self, reaction, lb=None, ub=None):
        """
        Modifies the allowed flux through the reaction(s).

        The method can also be used to modify multiple reactions at the same time.
        In that case each argument must be an iterable of equal length and
        matching positions corresponding to the same reaction.

        Parameters
        ----------
        reaction: iterable or `BasicReaction`
            A single `BasicReaction` instance or an iterable with multiple ones.
        lb: float (optional)
            A lower bound on the mass flux through this reaction. A default
            value can be set in the options.
        ub: float (optional)
            An upper bound on the mass flux through this reaction. A default
            value can be set in the options.
        """
        raise self._error

    def modify_compound_bounds(self, compound, lb=None, ub=None):
        """
        Modifies the bounds on compounds.

        Parameters
        ----------
        compound: iterable or `BasicCompound`
            A single `BasicCompound` instance or an iterable with multiple ones.
        lb: float (optional)
            A lower bound on the mass flux through this reaction. A default
            value can be set in the options.
        ub: float (optional)
            An upper bound on the mass flux through this reaction. A default
            value can be set in the options.
        """
        raise self._error

    def get_reaction_bounds(self, reaction=None):
        """
        Query the bounds on reactions currently in place.

        Parameters
        ----------
        reaction: iterable or `BasicReaction`
            A single `BasicReaction` instance or an iterable with multiple ones.

        Returns
        -------
        tuple or iterator:
            The lower and upper bound on a reaction or an iterator over pairs in
            the same order as the iterable provided. With no argument the
            returned iterator walks over triples of all present reactions with
            their respective lower and upper bounds.
        """
        raise self._error

    def add_compound_drain(self, compound, lb=None, ub=None):
        """
        Adds a drain for a certain compound or compounds to the model.

        Parameters
        ----------
        compound: iterable or `BasicCompound`
            A single `BasicCompound` instance or an iterable with multiple ones.
        lb: float (optional)
            A lower bound on the mass flux through this reaction. A default
            value can be set in the options.
        ub: float (optional)
            An upper bound on the mass flux through this reaction. A default
            value can be set in the options.
        """
        raise self._error

    def add_compound_source(self, compound, lb=None, ub=None):
        """
        Adds a source for a certain compound or compounds to the model.

        Parameters
        ----------
        compound: iterable or `BasicCompound`
            A single `BasicCompound` instance or an iterable with multiple ones.
        lb: float (optional)
            A lower bound on the mass flux through this reaction. A default
            value can be set in the options.
        ub: float (optional)
            An upper bound on the mass flux through this reaction. A default
            value can be set in the options.
        """
        raise self._error

    def get_reactions(self):
        """
        Returns
        -------
        iterator:
            Iterator over all reactions (excluding sources and drains).
        """
        raise self._error

    def get_compounds(self):
        """
        Returns
        -------
        iterator:
            Iterator over all compounds.
        """
        raise self._error

    def get_drains(self):
        """
        Returns
        -------
        iterator:
            Iterator over all drains.
        """
        raise self._error

    def get_sources(self):
        """
        Returns
        -------
        iterator:
            Iterator over all sources.
        """
        raise self._error

    def get_substrates_and_products(self, reaction, coefficients=False):
        """
        Parameters
        ----------
        reaction: iterable or `BasicReaction`
            A single `BasicReaction` instance or an iterable with multiple ones.
        coefficients: bool (optional)

        Returns
        -------
        iterator:
            Iterator over `BasicCompound` instances or pairs with coefficients.
            Returns an iterator over iterators of these if `reaction` is
            iterable.
        """
        raise self._error

    def get_objective_reaction(self, factor=False):
        """
        Parameters
        ----------
        factor: bool (optional)
            Causes the returned iterator to run over pairs of reaction and
            weight in the objective function.

        Returns
        -------
        iterator:
            Current reaction(s) that are used as objectives in LP.
        """
        raise self._error

    def set_objective_reaction(self, reaction, factor):
        """
        Sets a certain reaction as objective (for FBA).

        Parameters:
        -------
        reaction: iterable or `BasicReaction`
            A single `BasicReaction` instance or an iterable with multiple ones.
        factor: iterable or float
            Weight of the reaction in the objective function.
        """
        raise self._error

    def get_objective_value(self, threshold=None):
        """
        Parameters
        ----------
        threshold: float (optional)
            Value below which the objective value is considered to be zero. By
            default the model precision is used.

        Returns
        -------
        float:
            Flux of the set objective reaction(s).
        """
        raise self._error

    def fba(self, maximize=True):
        """
        Performs an optimization of the current objective(s) in the model.
        """
        raise self._error

    def parsimonious_fba(self):
        """
        Performs an optimization of the current objective(s) in the model
        followed by a minimisation of all other unnecessary fluxes.
        """
        raise self._error

    def get_flux_distribution(self, threshold=None):
        """
        Parameters
        ----------
        threshold: float (optional)
            Value below which a flux value is considered to be zero. By
            default the model precision is used.

        Returns
        -------
        iterator:
            Iterator over pairs of `BasicReaction`s and their flux.
        """
        raise self._error

    def set_medium(self, medium, lb=None, ub=None):
        """
        Modifies the allowed flux of the compound sources thereby specifying the
        growth medium composition.

        Parameters
        ----------
        medium: iterable
            Iterable of `BasicCompounds` that should in the growth medium.
        lb: iterable (optional)
            Lower bounds on the mass flux of the sources. A default
            value can be set in the options.
        ub: iterable (optional)
            Upper bounds on the mass flux of the sources. A default
            value can be set in the options.
        """
        raise self._error

    def is_fixed(self, reaction):
        """
        Tests whether a reaction's lower and upper bound are equal.

        Parameters
        ----------
        reaction: iterable or `BasicReaction`
            A single `BasicReaction` instance or an iterable with multiple ones.
        """
        raise self._error

    def export2lp(self, filename):
        """
        Export the current model to a flat text file in *.lp format.
        """
        raise self._error



#class FBAModel(object):
#    """
#    """
#
#    def __init__(self, model):
#        """
#        description
#
#        Parameters
#        ----------
#        name: type:
#            description
#
#        Returns
#        -------
#        type:
#            description
#
#        Notes
#        -----
#        text
#        """
#        self._model = model
#
#    def __copy__(self):
#        return FBAModel(copy.copy(self._model))
#
#    def __deepcopy__(self, memo=dict()):
#        return FBAModel(copy.deepcopy(self._model))
#
#    def copy(self):
#        return FBAModel(copy.deepcopy(self._model))
#
#    def add_reaction(self, name, coefficients, bounds=tuple()):
#        self._model.add_column(name, coefficients, bounds)
#
#    def knockout_reaction(self, reaction):
#        """
#        Constrains the allowed reaction flux to zero.
#
#        Parameters
#        ----------
#        reaction: iterable or str
#            A single name or iterable that contains the names of the reactions.
#        """
#        if hasattr(reaction, "__iter__"):
#            bounds = dict((rxn, (0.0, 0.0)) for rxn in reaction)
#        else:
#            bounds = {reaction: (0.0, 0.0)}
#        self._model.modify_column_bounds(bounds)
#
#    def delete_reaction_from_stoich(self, reaction):
#        """
#        Completely removes the reaction(s) from the model.
#
#        Parameters
#        ----------
#        reaction: iterable or str
#            A single name or iterable that contains the names of the reactions.
#        """
#        self._model.delete_column(reaction)
#
#    def knockout_compound(self, compound):
#        """
#        Constrains the reaction fluxes consuming or producing elements in
#        compound to zero.
#
#        Parameters
#        ----------
#        compound: iterable or str
#            A single name or iterable that contains the names of the compounds.
#        """
#        if hasattr(compound, "__iter__"):
#            bounds = dict((name, (0.0, 0.0)) for cmpd in compound for name in
#                    self._model.get_column_names(compound))
#        else:
#            bounds = dict((name, (0.0, 0.0)) for name in
#                    self._model.get_column_names(compound))
#        self._model.modify_column_bounds(metabolites)
#
#    def delete_compound_from_stoich(self, compound):
#        """
#        Completely removes the compound(s) from the model.
#
#        Parameters
#        ----------
#        reaction: iterable or str
#            A single name or iterable that contains the names of the compounds.
#        """
#        self._model.delete_row(compound)
#
#    def free_metabolites(self, compound):
#        """
#        Completely removes any bounds on compound(s).
#
#        Parameters
#        ----------
#        compound: iterable or str
#            A single name or iterable that contains the names of the compounds.
#        """
#        if hasattr(compound, "__iter__"):
#            bounds = dict((name, (-numpy.inf, numpy.inf)) for name in compound)
#        else:
#            bounds = {compound: (-numpy.inf, numpy.inf)}
#        self._model.modify_row_bounds(boundaries)
#
#    def modify_reaction_bounds(self, bounds):
#        """
#        Modifies the allowed flux on reactions.
#
#        Parameters
#        ----------
#        bounds : dict
#            Map of reaction names to pairs with their lower and upper bound.
#        """
#        self._model.modify_column_bounds(bounds)
#
#    def modify_compound_bounds(self, bounds):
#        """
#        Modifies the bounds on compounds.
#
#        Parameters
#        ----------
#        bounds : dict
#            Map of compound names to pairs with their lower and upper bound.
#        """
#        self._model.modify_row_bounds(self, bounds)
#
#    def get_reaction_bounds(self, reaction):
#        """
#        Query the bounds on reactions currently in place.
#
#        Parameters
#        ----------
#        reaction : str or iterable
#            A reaction identifier or an iterable of reaction identifiers.
#
#        Returns
#        -------
#        tuple or iterator over tuples:
#            The lower and upper bound on a reaction or an iterator over pairs in
#            the same order as the iterable provided.
#        """
#        if hasattr(reaction, "__iter__"):
#            return (self._model.get_column_bounds(rxn) for rxn in reaction)
#        else:
#            return self._model.get_column_bounds(reaction)
#
#    def add_compound_drain(self, compound, suffix="_Drain", bounds=tuple()):
#        """
#        Adds an export channel for each compound.
#
#        Parameters
#        ----------
#        compound: iterable or str
#            A single name or iterable that contains the names of the compounds
#            for which export from the model is added.
#        suffix: str (optional)
#            Identifier added to the name of drains to distinguish them from
#            normal reactions.
#        bounds: tuple (optional)
#            A pair of lower and upper bound on the flux for all drains added.
#        """
#        if hasattr(compound, "__iter__"):
#            columns = ((name + suffix, {name: -1.0}, bounds) for name in compound)
#            self._model.add_columns(columns)
#        else:
#            self._model.add_column(compound + suffix, {compound: -1.0}, bounds)
#
#    def add_compound_transport(self, compound, suffix = "_Transp", bounds=tuple()):
#        """
#        Adds an export channel for each compound.
#
#        Parameters
#        ----------
#        compound : iterable or str
#            A single name or iterable that contains the names of the compounds
#            for which export from the model is added.
#        suffix: str (optional)
#            Identifier added to the name of transporters to distinguish them from
#            normal reactions.
#        bounds: tuple (optional)
#            A pair of lower and upper bound on the flux for all transporters added.
#        """
#        if hasattr(compound, "__iter__"):
#            columns = ((name + suffix, {name: 1.0}, bounds) for name in compound)
#            self._model.add_columns(columns)
#        else:
#            self._model.add_column(compound + suffix, {compound: 1.0}, bounds)
#
#    def get_reactions(self, drain="_Drain", transp="_Transp"):
#        """
#        Returns
#        -------
#        iterator:
#            Iterator over all names of reactions (excluding transport and drain
#            reactions).
#        """
#        return (rxn for rxn in self._model.get_column_names() if not
#                (rxn.endswith(drain) or rxn.endswith(transp)))
#
#    def get_compounds(self):
#        """
#        Returns
#        -------
#        iterator:
#            Iterator over all names of compounds.
#        """
#        return self._model.get_row_names()
#
#    def get_drains(self, drain="_Drain"):
#        """
#        Gets the list of all transporters.
#
#        Returns
#        -------
#        iterator:
#            Iterator over all names of drains.
#        """
#        return (rxn for rxn in self._model.get_column_names() if
#                rxn.endswith(drain))
#
#    def get_transporters(self, transp="_Transp"):
#        """
#        Gets the list of all transporters.
#
#        Returns
#        -------
#        iterator:
#            Iterator over all names of transporters.
#        """
#        return (rxn for rxn in self._model.get_column_names() if
#                rxn.endswith(transp))
#
#    def get_substrates_and_products(self, reaction, coefficients=False):
#        """
#        Parameters
#        ----------
#        coefficients: bool (optional)
#            Causes the returned lists to include pairs with compound names and
#            coefficients.
#
#        Returns
#        -------
#        tuple:
#            Pair of the names of the substrates and the names of the products
#            of the reaction as lists.
#        """
#
#        def compound_sorting():
#            if factor > 0.0:
#                products.append(cmpd)
#            elif factor < 0.0:
#                substrates.append(cmpd)
#
#        def compound_sorting_with_factors():
#            if factor > 0.0:
#                products.append((cmpd, factor))
#            elif factor < 0.0:
#                substrates.append((cmpd, factor))
#
#        substrates=list()
#        products=list()
#        if coefficients:
#            sort = compound_sorting_with_factors
#        else:
#            sort= compound_sorting
#        for (cmpd, factor) in self._model.get_row_names(reaction, True):
#            sort()
#        return (substrates,products)
#
#    def get_objective_reaction(self, coefficients=False):
#        """
#        Parameters
#        ----------
#        coefficients: bool (optional)
#            Causes the returned iterator to run over pairs of reaction name and
#            absolute weight in the objective function.
#
#        Returns
#        -------
#        iterator:
#            Current reaction(s) that are used as objectives in LP.
#        """
#        return self._model.get_objective(coefficients)
#
#    def set_objective_reaction(self, reaction):
#        """
#        Sets a certain reaction as objective (for maximization).
#
#        Parameters:
#        -------
#        reaction: dict
#            Map from reaction name(s) to their factor(s).
#        """
#        self._model.set_objective(reaction)
#
#    def get_objective_value(self, threshold=1E-6):
#        """
#        Parameters
#        ----------
#        threshold: float (optional)
#            Value below which the objective value is considered to be zero.
#
#        Returns
#        -------
#        float:
#            Flux of the set objective reaction(s).
#        """
#        growth = self._model.get_objective_value()
#        return growth if growth > threshold else 0.0
#
##    def set_objective_reaction_minimize_rest(self, reaction, factor):
##        """
##        Sets a certain reaction as objective (for maximization) and minimizes
##        the usage of all the other reactions. All the transporter reactions are
##        set to zero flux.
##
##        Parameters:
##        -------
##        reaction: iterable
##            description: iterable that contains the name of the reaction.
##        factor: iterable
##            description: iterable that contains the factor to which the fluxes
##            of the reactions that are not objective should be minimized
##        """
##        react_dict=dict()
##        react_list=self._model.get_column_names()
##        transp_list=self.get_transporters()
##        for r in react_list:
##            if r in transp_list:
##                react_dict[r]=(0,0)
##            else:
##                react_dict[r]=(factor,factor)
##        self.modify_reaction_bounds(react_dict)
##        self.set_reaction_objective(reaction)
#
#    def fba(self, maximize=True):
#        """
#        Performs an optimization of the current objective(s) in the model.
#        """
#        return self._model.optimize(maximize)
#
#    def get_flux_distribution(self):
#        return self._model.get_solution_vector()
#
#    def set_medium(self, medium, upper = 20.0, transp="_Transp"):
#        """
#        Applies the given medium to the model by setting the corresponding
#        transporter bounds of all components to 0, upper
#        """
#        # reset all current transporter boundaries as a safety measure
#        bounds = dict(itertools.izip(self.get_transporters(transp),
#                itertools.repeat((0.0, 0.0))))
#        self.modify_reaction_bounds(bounds)
#        bounds = dict(itertools.izip(medium, itertools.repeat((0.0, upper))))
#        self.modify_reaction_bounds(bounds)
##        self.medium = medium #we could also read this from the bounds, but this way is much easier
#
#    def is_fixed(self, reaction):
#        """
#        Tests whether a reaction's lower and upper bound are equal.
#        """
#        bounds = self.get_reaction_bounds(reaction)
#        return bounds[0] == bounds[1]
#
#    def export2lp(self, filename):
#        self._model.export2lp(filename)


def generate_random_medium(transporters, percentage_range=(5, 100),
        minimal=list(), transp="_Transp"):
    """
    Generates a completely random medium based on a percentage of activated
    transporters.

    Parameters:
    -------
    transporters:
        asdfd
    percentage_range: tuple
        A random percentage of transporters is considered for the random medium
        according to this range. The first of the pair must be smaller than or
        equal to the second.
    minimal: iterable
        Some always active transporters that form a minimal medium that is
        extended with random other components.
    transp: str
        The suffix for transporters in the model.
    """
    assert percentage_range[0] <= percentage_range[1]

    # only choose from non-minimal transporters
    # -> ensures constant active percentage by preventing overlap
    choices = [t for t in transporters if not t in minimal]

    # select a random percentage of active transporters
    active = random.sample(choices, int(numpy.ceil(len(choices) *
            random.uniform(*percentage_range) / 100.0)))

    # since bounds is a dictionary we do not care about duplicates
    for trns in minimal:
        active.append(trns)

    return active

def set_random_medium(model, default_bound=(20.0, 20.0),
        percentage_range=(5, 100), minimal=list(), transp="_Transp"):
    """
    Generates and sets a completely random medium based on a percentage
    of activated transporters.

    Parameters:
    -------
    model: FBAModel
        The model is modified directly.
    default_bound: tuple
        A default upper limit for all components of the medium is chosen from
        this range of floating point numbers. The first of the pair must be
        smaller than or equal to the second.
    percentage_range: tuple
        A random percentage of transporters is considered for the random medium
        according to this range. The first of the pair must be smaller than or
        equal to the second.
    minimal: iterable
        Some always active transporters that form a minimal medium that is
        extended with random other components.
    transp: str
        The suffix for transporters in the model.
    """
    assert default_bound[0] <= default_bound[1]

    medium = generate_random_medium(list(model.get_transporters()),
                                     percentage_range, minimal, transp)
    upper = random.uniform(*default_bound)

    return model.set_medium(medium, upper)

