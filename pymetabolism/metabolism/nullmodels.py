#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=====================
Metabolic Null Models
=====================

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-07-01
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    nullmodels.py
"""


import itertools
import numpy
import networkx as nx

from decimal import Decimal
from collections import defaultdict
from Queue import Queue


#def bfs_balance(network, coefficients, seed=None):
#    """
#    Build a consistent stoichiometric matrix from a directed bipartite graph
#    object.
#
#    Requires a distribution of stoichiometric coefficients that will be
#    randomly chosen to construct mass-balanced reactions. Note that the
#    stoichiometric coefficients can be biased by inserting certain numbers
#    multiple times.
#
#    Basically, this function performs a breadth-first search on the network
#    assigning stoichiometric coefficients to reactions whenever possible and
#    otherwise reconciles compound masses.
#
#    Parameters
#    ----------
#    network: MetabolicNetwork
#        A bipartite MetabolicNetwork with a desired architecture.
#    coefficients: iterable
#        An iterable of stoichiometric coefficients that may be weighted by
#        including certain coefficients multiple times.
#    seed: int (optional)
#        A specific seed for the random number generator for reproducible runs.
#    """
#    rand_int = numpy.random.random_integers
#    if seed:
#        numpy.random.seed(int(seed))
#    # a map that stores exact masses of metabolites, thus decimal type
#    mass = defaultdict(dec.Decimal) # zero initialised decimal values
#    # stoichiometric matrix
#    stoichiometry = numpy.zeros(shape=(len(network.compounds),
#        len(network.reactions)), dtype=int)
#    # breadth-first search containers
#    disc = Queue() # queue for discovered nodes
#    nbunch = list(network.compounds)
#    # initialise search
#    elem = nbunch[rand_int(0, len(nbunch) - 1)]
#    mass[elem] = Decimal(1)
#    disc.append(elem)
#    while True:
#        elem = disc.pop()
#        if mass[elem] > 0:
#            continue
#        for node in graph.successors_iter(elem):
#            disc.append(node)
#    # do some consistency checks for example all masses > 0 and
#    # each reaction must satisfy mass substrates = mass products

def lp_balance(network, coefficients, solver="gurobi", seed=None):
    """
    Build a consistent stoichiometric matrix from a directed bipartite graph
    object.

    Requires a distribution of stoichiometric coefficients that will be
    randomly chosen to construct mass-balanced reactions. Note that the
    stoichiometric coefficients can be biased by inserting certain numbers
    multiple times.

    Mass-balance is achieved by solving the following linear programming
    problem:
        Minimise the sum over all stoichiometric coefficients subject to the
        constraints given by the transpose of the stoichiometric matrix
        multiplied with the mass vector equals the zero vector.

    Parameters
    ----------
    network: MetabolicNetwork
        A bipartite MetabolicNetwork with a desired architecture.
    coefficients: iterable
        An iterable of stoichiometric coefficients that may be weighted by
        including certain coefficients multiple times.
    seed: int (optional)
        A specific seed for the random number generator for reproducible runs.
    """
    if seed:
        numpy.random.seed(int(seed))
    if seed:
        rnd.seed(seed)
    # a map that stores exact masses of metabolites, thus decimal type
    mass = dict()
    # map with metabolite matrix-indeces
    metb_idx = dict()
    for (i, metb) in enumerate(metbs):
        metb_idx[metb] = i
        mass[metb] = rnd.random()
#    # map with reaction rates
#    rates = dict()
    # map with reaction matrix-indeces
    rxn_idx = dict()
    for (i, rxn) in enumerate(rxns):
        rxn_idx[rxn] = i
#        rates[rxn] = rnd.random() # not required atm
    m = len(metbs)
    n = len(rxns)
    # problem: minimise sum over stoichiometric coefficients
    f = np.ones([n, m], dtype=float)
    # subject to T(S).m = 0
    a_eq = np.empty([n, m], dtype=float)
    for rxn in rxns:
        for metb in metbs:
            a_eq[rxn_idx[rxn], metb_idx[metb]] = mass[metb]
    b_eq = np.zeros(n, dtype=float)
    # where
    lb = np.zeros([n, m], dtype=float)
    for rxn in rxns:
        for metb in graph.predecessors_iter(rxn):
            lb[rxn_idx[rxn], metb_idx[metb]] = -100.
        for metb in graph.successors_iter(rxn):
            lb[rxn_idx[rxn], metb_idx[metb]] = 1.
    ub = np.zeros([n, m], dtype=float)
    for rxn in rxns:
        for metb in graph.predecessors_iter(rxn):
            lb[rxn_idx[rxn], metb_idx[metb]] = -1.
        for metb in graph.successors_iter(rxn):
            lb[rxn_idx[rxn], metb_idx[metb]] = 100.
    # solve
    p = oo.LP(f=f, A=None, Aeq=a_eq, b=None, beq=b_eq, lb=lb, ub=ub)
    #p.debug = 1
    result = p.solve('cvxopt_lp')
    print result.ff
    print result.xf

