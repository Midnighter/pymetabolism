#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
============================
Metabolic Network Generators
============================

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-07-01
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    generators.py
"""

import logging
import numpy
import numpy.random

from ..metabolism import metabolism as met
from ..network import networks as nets
from .. import miscellaneous as misc


logger = logging.getLogger(__name__)
logger.addHandler(misc.NullHandler())


def prune_network(network):
    """
    Removes stub reactions (in- or out-degree of 1 and 0 respectively) and
    assures that all other reactions consume and produce something.

    Parameters
    ----------
    network: MetabolicNetwork
        A MetabolicNetwork instance.
    """
    rand_int = numpy.random.random_integers
    num_rxns = len(network.reactions)
    num_cmpds = len(network.compounds)
    total = 0
    prune = list()
    for rxn in network.reactions:
        in_deg = network.in_degree(rxn)
        out_deg = network.out_degree(rxn)
        if in_deg == 0:
            if out_deg <= 1:
                prune.append(rxn)
#                network.remove_node(rxn)
#                logger.debug("removed reaction %s", str(rxn))
            else:
                targets = network.successors(rxn)
                flips = rand_int(1, len(targets) - 1)
                while (flips > 0):
                    target = targets[rand_int(0, len(targets) - 1)]
                    factor = network[rxn][target]["coefficient"]
                    network.remove_edge(rxn, target)
                    network.add_edge(target, rxn, coefficient=factor)
                    logger.debug("flipped direction of link %s -> %s",
                            str(rxn), str(target))
                    targets.remove(target)
                    flips -= 1
                    total += 1
        elif out_deg == 0:
            if in_deg <= 1:
                prune.append(rxn)
#                network.remove_node(rxn)
#                logger.debug("removed reaction %s", str(rxn))
            else:
                targets = network.predecessors(rxn)
                flips = rand_int(1, len(targets) - 1)
                while (flips > 0):
                    target = targets[rand_int(0, len(targets) - 1)]
                    factor = network[target][rxn]["coefficient"]
                    network.remove_edge(target, rxn)
                    network.add_edge(rxn, target, coefficient=factor)
                    logger.debug("flipped direction of link %s -> %s",
                            str(rxn), str(target))
                    targets.remove(target)
                    flips -= 1
                    total += 1
    for rxn in prune:
        network.remove_node(rxn)
        logger.debug("removed reaction %s", str(rxn))
    for cmpd in network.compounds:
        if network.degree(cmpd) == 0:
            network.remove_node(cmpd)
            logger.debug("removed compound %s", str(cmpd))
    logger.info("%d reactions and %d compounds removed",
            (num_rxns - len(network.reactions)),
            (num_cmpds - len(network.compounds)))
    logger.info("direction of %d links reversed", total)


def random_p_mn(num_compounds, num_reactions, num_reversible, p, seed=None):
    """
    Creates a bipartite graph that models metabolism according to the principles
    of an Erdos-Renyi-like random graph.

    Parameters
    ----------
    num_compounds: int
        The number of compounds (approximately) that should be in the network.
    num_reactions: int
        The number of reactions (approximately) that should be in the network.
    num_reversible: int
        The number of reactions that are reversible.
    p: float
        The probability that a link between a compound and a reaction exists.
    seed: int (optional)
        A specific seed for the random number generator for reproducible runs.

    Returns
    -------
    pymetabolism.network.networks.MetabolicNetwork

    Notes
    -----
    The numbers of compounds and reactions may be slightly less than desired
    because isolated nodes are removed from the network.
    """
    # setup
    rand_float = numpy.random.random_sample
    rand_int = numpy.random.random_integers
    if seed:
        numpy.random.seed(int(seed))
    num_compounds = int(num_compounds)
    num_reactions = int(num_reactions)
    num_reversible = int(num_reversible)
    p = float(p)
    options = misc.OptionsManager()
    network = nets.MetabolicNetwork()
    # add compounds
    for i in range(num_compounds):
        network.add_node(met.BasicCompound("%s%d" % (options.compound_prefix, i)))
    # choose a number of reactions as reversible
    reversibles = set()
    while len(reversibles) < num_reversible:
        reversibles.add(rand_int(0, num_reactions - 1))
    for i in range(num_reactions):
        if i in reversibles:
            network.add_node(met.BasicReaction(
                    "%s%d" % (options.reaction_prefix, i), reversible=True))
        else:
            network.add_node(met.BasicReaction(
                "%s%d" % (options.reaction_prefix, i)))
    for src in network.compounds:
        for tar in network.reactions:
            if rand_float() < p:
                network.add_edge(src, tar, coefficient=0)
                logger.debug("added link %s -> %s.", str(src), str(tar))
            # a conditional case here (elif not if) because we cannot determine
            # substrates and products from bidirectional edges
            elif rand_float() < p:
                network.add_edge(tar, src, coefficient=0)
                logger.debug("added link %s -> %s.", str(tar), str(src))
    prune_network(network)
    return network

def random_scale_free_mn(num_compounds, num_reactions, num_reversible,
        num_rxn_tar, num_cmpd_tar, seed=None):
    """
    Uses a Barabasi-Alberts-like preferential attachment algorithm. Adopted from
    the networkx implementation.

    Parameters
    ----------
    num_compounds: int
        The number of compounds that should be in the network.
    num_reactions: int
        The number of reactions that should be in the network.
    num_reversible: int
        The number of reactions that are reversible.
    num_rxn_tar: int
        How many compounds a new reaction node links to.
    num_cmpd_tar: int
        How many reactions a new compound node links to.
    seed: int (optional)
        A specific seed for the random number generator for reproducible runs.
    """
    # setup
    rand_int = numpy.random.random_integers
    rand_float = numpy.random.random_sample
    if seed:
        numpy.random.seed(int(seed))
    num_compounds = int(num_compounds)
    num_reactions = int(num_reactions)
    num_reversible = int(num_reversible)
    num_rxn_tar = int(num_rxn_tar)
    num_cmpd_tar = int(num_cmpd_tar)
    options = misc.OptionsManager()
    network = nets.MetabolicNetwork()
    # target nodes for reactions
    rxn_targets = []
    for i in range(num_rxn_tar):
        comp = met.BasicCompound("%s%d" % (options.compound_prefix, i))
        network.add_node(comp)
        rxn_targets.append(comp)
    # target nodes for compounds
    cmpd_targets = []
    # biased lists for preferential attachment
    repeated_cmpds = []
    repeated_rxns = []
    # choose a number of reactions as reversible
    reversibles = set()
    while len(reversibles) < num_reversible:
        reversibles.add(rand_int(0, num_reactions - 1))
    for i in range(num_cmpd_tar):
        if i in reversibles:
            rxn = met.BasicReaction("%s%d" % (options.reaction_prefix, i),
                    reversible=True)
        else:
            rxn = met.BasicReaction("%s%d" % (options.reaction_prefix, i))
        network.add_node(rxn)
        cmpd_targets.append(rxn)
        for cmpd in rxn_targets:
            if rand_float() < 0.5:
                network.add_edge(rxn, cmpd, coefficient=0)
                logger.debug("added link %s -> %s", str(rxn), str(cmpd))
            else:
                network.add_edge(cmpd, rxn, coefficient=0)
                logger.debug("added link %s -> %s", str(cmpd), str(rxn))
        repeated_cmpds.extend(rxn_targets)
        repeated_rxns.extend([rxn] * num_rxn_tar)
#    logger.debug("Targets for compounds: %s", comp_targets)
    # current vertices being added
    current_rxn = num_cmpd_tar
    current_cmpd = num_rxn_tar
    while (current_cmpd < num_compounds or current_rxn < num_reactions):
        if current_cmpd < num_compounds:
            source = met.BasicCompound("%s%d" % (options.compound_prefix,
                    current_cmpd))
            network.add_node(source)
            for rxn in cmpd_targets:
                if rand_float() < 0.5:
                    network.add_edge(source, rxn, coefficient=0)
                    logger.debug("added link %s -> %s", str(source), str(rxn))
                else:
                    network.add_edge(rxn, source, coefficient=0)
                    logger.debug("added link %s -> %s", str(rxn), str(source))
            repeated_rxns.extend(cmpd_targets)
            repeated_cmpds.extend([source] * num_cmpd_tar)
            cmpd_targets = set()
            while len(cmpd_targets) < num_cmpd_tar:
                rxn = repeated_rxns[rand_int(0, len(repeated_rxns) - 1)]
                cmpd_targets.add(rxn)
            current_cmpd += 1
        if current_rxn < num_reactions:
            if current_rxn in reversibles:
                source = met.BasicReaction("%s%d" % (options.reaction_prefix,
                    current_rxn), reversible=True)
            else:
                source = met.BasicReaction("%s%d" % (options.reaction_prefix,
                    current_rxn))
            network.add_node(source)
            for cmpd in rxn_targets:
                if rand_float() < 0.5:
                    network.add_edge(source, cmpd, coefficient=0)
                    logger.debug("added link %s -> %s", str(source), str(cmpd))
                else:
                    network.add_edge(cmpd, source, coefficient=0)
                    logger.debug("added link %s -> %s", str(cmpd), str(source))
            repeated_cmpds.extend(rxn_targets)
            repeated_rxns.extend([source] * num_rxn_tar)
            rxn_targets = set()
            while len(rxn_targets) < num_rxn_tar:
                cmpd = repeated_cmpds[rand_int(0, len(repeated_cmpds) - 1)]
                rxn_targets.add(cmpd)
            current_rxn += 1
    prune_network(network)
    return network

