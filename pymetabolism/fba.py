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

    def initialise(self, name=""):
        """
        Initialise a new and empty model.

        Parameters
        ----------
        name: str (optional)
            Identify the model by a name.
        """
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
        self.modify_reaction_bounds(reaction, lb=0.0, ub=0.0)

    def free_reaction(self, reaction):
        """
        Completely removes any bounds on the reaction(s).

        Parameters
        ----------
        reaction: iterable or `BasicReaction`
            A single `BasicReaction` instance or an iterable with multiple ones.
        """
        self.modify_reaction_bounds(reaction, lb=-numpy.inf, ub=numpy.inf)

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

    def modify_reaction_coefficients(self, reaction, coefficients):
        """
        Modify the coefficients of compounds in the reaction.

        This can be used to introduce an existing compound into a reaction that
        it previously did not take part in.

        Parameters
        ----------
        reaction: `BasicReaction` or iterable
            An instance of `BasicReaction`.
        coefficients: iterable
            Pairs of `BasicCompound` instances and their
            stoichiometric coefficient.
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

    def iter_reaction_bounds(self, reaction=None):
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

    def iter_reactions(self, compound=None, coefficients=False):
        """
        Parameters
        ----------
        compound: `BasicCompound` (optional)
            A single `BasicCompound` instance.
        coefficients: bool (optional)
            Return also the stoichiometric coefficients of the given compound in
            the reactions it participates in.

        Returns
        -------
        iterator:
            Iterator over all reactions (excluding sources and drains) that
            `compound` is involved in.
        """
        raise self._error

    def iter_compounds(self):
        """
        Returns
        -------
        iterator:
            Iterator over all compounds.
        """
        raise self._error

    def iter_drains(self):
        """
        Returns
        -------
        iterator:
            Iterator over all drains.
        """
        raise self._error

    def iter_sources(self):
        """
        Returns
        -------
        iterator:
            Iterator over all sources.
        """
        raise self._error

    def iter_substrates_and_products(self, reaction, coefficients=False):
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

    def iter_objective_reaction(self, factor=False):
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

    def iter_flux(self, reaction=None, threshold=None):
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

    def set_medium(self, compound, lb=None, ub=None):
        """
        Modifies the allowed flux of the compound sources thereby specifying the
        growth medium composition.

        Parameters
        ----------
        compound: iterable or `BasicCompounds`
            Iterable of `BasicCompounds` that should be in the growth medium.
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


#def generate_random_medium(transporters, percentage_range=(5, 100),
#        minimal=list(), transp="_Transp"):
#    """
#    Generates a completely random medium based on a percentage of activated
#    transporters.
#
#    Parameters:
#    -------
#    transporters:
#        asdfd
#    percentage_range: tuple
#        A random percentage of transporters is considered for the random medium
#        according to this range. The first of the pair must be smaller than or
#        equal to the second.
#    minimal: iterable
#        Some always active transporters that form a minimal medium that is
#        extended with random other components.
#    transp: str
#        The suffix for transporters in the model.
#    """
#    assert percentage_range[0] <= percentage_range[1]
#
#    # only choose from non-minimal transporters
#    # -> ensures constant active percentage by preventing overlap
#    choices = [t for t in transporters if not t in minimal]
#
#    # select a random percentage of active transporters
#    active = random.sample(choices, int(numpy.ceil(len(choices) *
#            random.uniform(*percentage_range) / 100.0)))
#
#    # since bounds is a dictionary we do not care about duplicates
#    for trns in minimal:
#        active.append(trns)
#
#    return active
#
#def set_random_medium(model, default_bound=(20.0, 20.0),
#        percentage_range=(5, 100), minimal=list(), transp="_Transp"):
#    """
#    Generates and sets a completely random medium based on a percentage
#    of activated transporters.
#
#    Parameters:
#    -------
#    model: FBAModel
#        The model is modified directly.
#    default_bound: tuple
#        A default upper limit for all components of the medium is chosen from
#        this range of floating point numbers. The first of the pair must be
#        smaller than or equal to the second.
#    percentage_range: tuple
#        A random percentage of transporters is considered for the random medium
#        according to this range. The first of the pair must be smaller than or
#        equal to the second.
#    minimal: iterable
#        Some always active transporters that form a minimal medium that is
#        extended with random other components.
#    transp: str
#        The suffix for transporters in the model.
#    """
#    assert default_bound[0] <= default_bound[1]
#
#    medium = generate_random_medium(list(model.get_transporters()),
#                                     percentage_range, minimal, transp)
#    upper = random.uniform(*default_bound)
#
#    return model.set_medium(medium, upper)

