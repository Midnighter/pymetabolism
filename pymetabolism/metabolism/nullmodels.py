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


import logging
import itertools
import numpy

from ..errors import PyMetabolismError
from .. import miscellaneous as misc


logger = logging.getLogger(__name__)
logger.addHandler(misc.NullHandler())


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

#def lp_balance(network, coefficients, solver="gurobi", seed=None):
#    """
#    Build a consistent stoichiometric matrix from a directed bipartite graph
#    object.
#
#    Requires a distribution of stoichiometric coefficients that will be
#    randomly chosen to construct mass-balanced reactions. Note that the
#    stoichiometric coefficients can be biased by inserting certain numbers
#    multiple times.
#
#    Mass-balance is achieved by solving the following linear programming
#    problem:
#        Minimise the sum over all stoichiometric coefficients subject to the
#        constraints given by the transpose of the stoichiometric matrix
#        multiplied with the mass vector equals the zero vector.
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
#    if seed:
#        numpy.random.seed(int(seed))
#    if seed:
#        rnd.seed(seed)
#    # a map that stores exact masses of metabolites, thus decimal type
#    mass = dict()
#    # map with metabolite matrix-indeces
#    metb_idx = dict()
#    for (i, metb) in enumerate(metbs):
#        metb_idx[metb] = i
#        mass[metb] = rnd.random()
##    # map with reaction rates
##    rates = dict()
#    # map with reaction matrix-indeces
#    rxn_idx = dict()
#    for (i, rxn) in enumerate(rxns):
#        rxn_idx[rxn] = i
##        rates[rxn] = rnd.random() # not required atm
#    m = len(metbs)
#    n = len(rxns)
#    # problem: minimise sum over stoichiometric coefficients
#    f = np.ones([n, m], dtype=float)
#    # subject to T(S).m = 0
#    a_eq = np.empty([n, m], dtype=float)
#    for rxn in rxns:
#        for metb in metbs:
#            a_eq[rxn_idx[rxn], metb_idx[metb]] = mass[metb]
#    b_eq = np.zeros(n, dtype=float)
#    # where
#    lb = np.zeros([n, m], dtype=float)
#    for rxn in rxns:
#        for metb in graph.predecessors_iter(rxn):
#            lb[rxn_idx[rxn], metb_idx[metb]] = -100.
#        for metb in graph.successors_iter(rxn):
#            lb[rxn_idx[rxn], metb_idx[metb]] = 1.
#    ub = np.zeros([n, m], dtype=float)
#    for rxn in rxns:
#        for metb in graph.predecessors_iter(rxn):
#            lb[rxn_idx[rxn], metb_idx[metb]] = -1.
#        for metb in graph.successors_iter(rxn):
#            lb[rxn_idx[rxn], metb_idx[metb]] = 100.
#    # solve
#    p = oo.LP(f=f, A=None, Aeq=a_eq, b=None, beq=b_eq, lb=lb, ub=ub)
#    #p.debug = 1
#    result = p.solve('cvxopt_lp')
#    print result.ff
#    print result.xf

def make_consistent_stoichiometry(network, coefficients, masses=None):
    """
    Based on a given network architecture this function attempts to generate a
    consistent stoichiometry that obeys mass conservation laws.

    Parameters
    ----------
    network: MetabolicNetwork
    coefficients: indexable
    masses: indexable
    """
    rand_int = numpy.random.random_integers

    def balance_reaction_by_mass(reaction):
        """
        Balance a single reaction by adjusting the stoichiometric coefficients in a
        way that leads to mass conservation.
        """
        # modify the coefficients for the current reaction
        temp_coeff = factors.copy()
        # change the boundaries so only involved compounds are affected
        temp_bounds = bounds.copy()

        # substrates
        msg = list()
        index = 1
        for cmpd in network.pred[reaction]:
            temp_coeff[cmpd.name] = -mass_vector[cmpd]
            temp_bounds[cmpd.name] = (0.1, upper)
            msg.append("-%s S%d" % (str(mass_vector[cmpd]), index))
            index += 1
        #products
        for cmpd in network.succ[reaction]:
            temp_coeff[cmpd.name] = mass_vector[cmpd]
            temp_bounds[cmpd.name] = (0.1, upper)
            msg.append("+%s S%d" % (str(mass_vector[cmpd]), index))
            index += 1
        msg.append("= 0")
        logger.debug("%s:", reaction.name)
        logger.debug("%s", " ".join(msg))

        model.modify_row_coefficients("reaction", temp_coeff)
        model.modify_column_bounds(temp_bounds)

        model.optimize(maximize=False)
        try:
            coeffs = dict(model.get_solution_vector())
        except PyMetabolismError as err:
            logger.debug(str(err))
            raise PyMetabolismError("Reaction '%s' cannot be balanced with the"\
                " given mass vector." % reaction.name)

        for cmpd in network.pred[reaction]:
            network[cmpd][reaction]["coefficient"] = coeffs[cmpd.name]
        for cmpd in network.succ[reaction]:
            network[reaction][cmpd]["coefficient"] = coeffs[cmpd.name]

    options = misc.OptionsManager()
    if not masses:
        # the default masses for compounds 
        masses = numpy.array([pair[1] for pair in\
            network.degree_iter(network.compounds)])
    # generate mass vector for compounds
    mass_vector = dict(itertools.izip(network.compounds,
            masses[rand_int(low=0, high=len(masses) - 1,
                size=len(network.compounds))]))
    # prepare a single LP model for all reactions
    model = options.get_lp_model()
    model.add_columns(((cmpd.name, {}, (0.0, max(coefficients))) for cmpd in\
            network.compounds))
    factors = dict(itertools.izip(model.get_column_names(), itertools.repeat(0.0)))
    model.add_row("reaction", factors)
    # test different objective functions:
    # - only zeros leads to fast solutions (no objective function) that are
    #   close to the upper boundary
    # - all equal leads to slower solutions that are mostly one with a few
    #   variations to balance the reactions
    # - an objective function with entries of 1 / factor leads to very long
    #   solution times but much more varied entries
    # - best approach so far: no objective function, but use a starting point
    #   by randomly picking coefficients
    model.set_objective(dict(itertools.izip(model.get_column_names(),
            itertools.repeat(1.0))))
#    model._make_integer(model.get_column_names())
    bounds = dict(itertools.izip(model.get_column_names(),
        itertools.repeat((0.0, 0.0))))
    upper = max(coefficients)

    total = float(len(network.reactions))
    for (i, rxn) in enumerate(network.reactions):
        balance_reaction_by_mass(rxn)
        logger.info("%.2f %% complete.", float(i + 1) / total * 100.)


