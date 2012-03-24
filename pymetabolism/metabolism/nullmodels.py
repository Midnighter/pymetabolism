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

from operator import itemgetter
from ..errors import PyMetabolismError
from .. import miscellaneous as misc
from ..fba import FBAModel


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

def make_consistent_stoichiometry(network, coefficients, mass_vector=None):
    """
    Based on a given network architecture this function attempts to generate a
    consistent stoichiometry that obeys mass conservation laws.

    Parameters
    ----------
    network: MetabolicNetwork
    coefficients: indexable
    mass_vector: dict
    """

    def balance_reaction_by_mass(reaction):
        """
        Balance a single reaction by adjusting the stoichiometric coefficients in a
        way that leads to mass conservation.
        """
        model.modify_reaction_bounds(network.compounds, lb=0.0)
        # abuse free_compound to reset all coefficients
        model.free_compound("reaction")

        compounds = [cmpd for cmpd in itertools.chain(network.pred[reaction],
                network.succ[reaction])]
        # modify the coefficients for the current reaction
        temp_coeff = list()

        # substrates
        msg = list()
        for cmpd in network.pred[reaction]:
            temp_coeff.append((cmpd, -mass_vector[cmpd]))
            msg.append("- %.2f %s" % (mass_vector[cmpd], str(cmpd)))
        #products
        for cmpd in network.succ[reaction]:
            temp_coeff.append((cmpd, mass_vector[cmpd]))
            msg.append("+ %.2f %s" % (mass_vector[cmpd], str(cmpd)))
        msg.append("= 0")
        logger.debug("%s:", reaction.name)
        logger.debug("%s", " ".join(msg))

        model.set_objective_reaction(compounds, 1.0)
        model.modify_reaction_bounds(compounds, lb=1.0)
        model.modify_compound_coefficients("reaction", temp_coeff)

        model.fba(maximize=False)
        try:
            coeffs = dict(model.iter_flux())
        except PyMetabolismError as err:
            logger.debug(str(err))
            raise PyMetabolismError("Reaction '%s' cannot be balanced with the"\
                " given mass vector.", str(reaction))

        for cmpd in network.pred[reaction]:
            network[cmpd][reaction]["coefficient"] = coeffs[cmpd]
        for cmpd in network.succ[reaction]:
            network[reaction][cmpd]["coefficient"] = coeffs[cmpd]

    if not mass_vector:
        # the default masses for compounds:
        # * compounds sorted by degree have a mass of degree equal to the other
        #   end of that sorted list (similar to inverse of degree but integers)
        compound = itemgetter(0)
        degree = itemgetter(1)
        masses = [pair for pair in network.degree_iter(network.compounds)]
        masses.sort(key=degree)
        end = len(masses) - 1
        # generate mass vector for compounds
        mass_vector = dict((compound(pair), degree(masses[end - i]))\
                for (i, pair) in enumerate(masses))

    # prepare a single LP model for all reactions
    model = FBAModel("mass balance")
    for cmpd in network.compounds:
        cmpd.reversible = False
    model.add_reaction(network.compounds, [("reaction", 0.0)], lb=0.0, ub=max(coefficients))
    # test different objective functions:
    # * only zeros leads to fast solutions (no objective function) that are
    #   close to the upper boundary, this can be ameliorated by picking starting
    #   points for the variables from the given distribution of coefficients
    # * all equal leads to slower solutions that are mostly one with a few
    #   variations to balance the reactions
    # * an objective function with entries of 1 / factor leads to very long
    #   solution times but much more varied coefficients

    # integer variables because traditionally biochemical reactions are
    # multiples of chemical groups, i.e., that's what will be in the
    # coefficients
    model._make_integer(network.compounds)

    total = float(len(network.reactions))
    for (i, rxn) in enumerate(network.reactions):
        balance_reaction_by_mass(rxn)
        logger.info("%.2f %% complete.", float(i + 1) / total * 100.)


