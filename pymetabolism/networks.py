#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=================================
Metabolic Network Representations
=================================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
    Nikolaus Sonnenschein
:Date:
    2011-04-13
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    network.py
"""


import copy
import networkx as nx
import pymetabolism.metabolism as pymet

from pymetabolism.errors import Error


class CompoundCentricNetwork(nx.MultiDiGraph):
    """
    """

    def __init__(self, name=""):
        """
        """
        nx.MultiDiGraph.__init__(self, name=name)

    def count_subgraphs(self):
        wrapper = mw.MfinderWrapper(self.to_directed())
        wrapper.count_subgraphs()
        self.graph["mtf_counts"] = wrapper.mtf_counts

    def to_directed(self):
        """
        Return a copy with no multiple edges and no attributes.
        """
        return nx.DiGraph(self.edges_iter())


class ReactionCentricNetwork(nx.MultiDiGraph):
    """
    """

    def __init__(self, name=""):
        """
        """
        nx.MultiDiGraph.__init__(self, name=name)

    def count_subgraphs(self):
        wrapper = mw.MfinderWrapper(self.to_directed())
        wrapper.count_subgraphs()
        self.graph["mtf_counts"] = wrapper.mtf_counts

    def to_directed(self):
        """
        Return a copy with no multiple edges and no attributes.
        """
        return nx.DiGraph(self.edges_iter())


class MetabolicNetwork(nx.DiGraph):
    """
    """

    def __init__(self, name=""):
        """
        """
        nx.DiGraph.__init__(self, name=name)
        self.reactions = set()
        self.compounds = set()

    def add_edge(self, u, v, **kw_args):
        """
        """
        if isinstance(u, pymet.BasicReaction):
            self.reactions.add(u)
        elif isinstance(u, pymet.BasicCompound):
            self.compounds.add(u)
        else:
            raise Error("unidentified metabolic component")
        if isinstance(v, pymet.BasicReaction):
            self.reactions.add(v)
        elif isinstance(v, pymet.BasicCompound):
            self.compounds.add(v)
        else:
            raise Error("unidentified metabolic component")
        nx.DiGraph.add_edge(self, u, v, **kw_args)

    def add_node(self, n, **kw_args):
        """
        """
        if isinstance(n, pymet.BasicReaction):
            self.reactions.add(n, **kw_args)
        elif isinstance(n, pymet.BasicCompound):
            self.compounds.add(n, **kw_args)
        else:
            raise Error("unidentified metabolic component")
        nx.DiGraph.add_node(self, n)

    def introduce_bidirectional(self):

        def add_rev(u, v):
            template.add_edge(u, v)
            template.add_edge(v, u)

        def add_single(u, v):
            template.add_edge(u, v)

        template = MetabolicNetwork(self.name)
        for cmpd in self.compounds:
            template.add_node(cmpd)
        for rxn in self.reactions:
            template.add_node(rxn)
        # introduce bidirectional edges
        for rxn in self.reactions:
            if rxn.reversible:
                new_edge = add_rev
            else:
                new_edge = add_single
            for cmpd in self.predecessors_iter(rxn):
                new_edge(cmpd, rxn)
            for cmpd in self.successors_iter(rxn):
                new_edge(rxn, cmpd)
        return template

    def count_subgraphs(self):
        raise NotImplementedError("motifs are not defined for bipartite graphs")

    def read_edgelist(self, path, compound_prefix="M", reaction_prefix="R_",
            reversible_suffix="_Rev", delimiter=None, comments="#"):
        """
        """

        def build_node(name):
            if name.startswith(compound_prefix):
                compound = pymet.BasicCompound(name[len(compound_prefix):])
#                self.compounds.add(compound)
                return compound
            elif name.startswith(reaction_prefix):
                if name.endswith(reversible_suffix):
                    reaction = pymet.BasicReaction(name[len(reaction_prefix):
                            -len(reversible_suffix)], reversible=True)
                    reaction.reversible = True
                else:
                    reaction = pymet.BasicReaction(name[len(reaction_prefix):])
#                self.reactions.add(reaction)
                return reaction
            else:
                raise TypeError("unrecognised metabolic object")

        with open(path, "r") as file_handle:
            lines = [line.strip() for line in file_handle]
        for line in lines:
            if line.startswith(comments) or line == "":
                continue
            tmp = line.split(delimiter)
            u = build_node(tmp[0])
            if isinstance(u, pymet.BasicReaction) and tmp[0].endswith(reversible_suffix):
                continue
            v = build_node(tmp[1])
            if isinstance(v, pymet.BasicReaction) and tmp[1].endswith(reversible_suffix):
                continue
            self.add_edge(u, v)

    def write_edgelist(self, path, distinct=True, compound_prefix="M_",
            reaction_prefix="R_", reversible_suffix="_Rev", delimiter="\t",
            comments="#"):
        """
        """
        lines = list()
        for rxn in self.reactions:
            rxn_name = reaction_prefix + rxn.name
            if rxn.reversible:
                if distinct:
                    rev_name = reaction_prefix + rxn.name + reversible_suffix
                else:
                    rev_name = rxn_name
                for cmpd in network.successors_iter(rxn):
                    lines.append("%s%s%s\n" % (rxn_name, delimiter, compound_prefix
                            + cmpd.name))
                    lines.append("%s%s%s\n" % (compound_prefix + cmpd.name,
                            delimiter, rev_name))
                for cmpd in network.predecessors_iter(rxn):
                    lines.append("%s%s%s\n" % (compound_prefix + cmpd.name,
                            delimiter, rxn_name))
                    lines.append("%s%s%s\n" % (rev_name, delimiter, compound_prefix
                            + cmpd.name))
            else:
                for cmpd in network.successors_iter(rxn):
                    lines.append("%s%s%s\n" % (rxn_name, delimiter, compound_prefix
                            + cmpd.name))
                for cmpd in network.predecessors_iter(rxn):
                    lines.append("%s%s%s\n" % (compound_prefix + cmpd.name,
                            delimiter, rxn_name))
        with open(path, "w") as file_handle:
            file_handle.writelines(lines)

    def to_compound_centric(self):
        """
        """

        def add_bi(u, v):
            network.add_edge(u, v)
            network.add_edge(v, u)

        network = CompoundCentricNetwork("compound_centric_" + self.name)
        # project to unipartite network with only compound nodes
        for cmpd in self.compounds:
            network.add_node(cmpd)
        for rxn in self.reactions:
            if rxn.reversible:
                # add a bidirectional link
                add_link = add_bi
            else:
                # add a unidirectional link
                add_link = network.add_edge
            for pred in self.predecessors_iter(rxn):
                for succ in self.successors_iter(rxn):
                    add_link(pred, succ)
        return network

    def to_reaction_centric(self):
        """
        """
        network = ReactionCentricNetwork("reaction_centric_" + self.name)
        # project to unipartite network with only reaction nodes
        for rxn in self.reactions:
            network.add_node(rxn)
        for cmpd in self.compounds:
            predecessors = self.predecessors(cmpd)
            rev_pred = [rxn for rxn in predecessors if rxn.reversible]
            successors = self.successors(cmpd)
            rev_succ = [rxn for rxn in successors if rxn.reversible]
            for pred in predecessors:
                for succ in successors:
                    network.add_edge(pred, succ)
            # add links due to reversibility
                for rxn in rev_pred:
                    network.add_edge(pred, rxn)
            for rxn in rev_succ:
                for succ in successors:
                    network.add_edge(rxn, succ)
                for pred in rev_pred:
                    network.add_edge(rxn, pred)
        # we added a lot of self-links in the process, I felt removing them
        # later was more efficient than working with set differences all the
        # time
        network.remove_edges_from(network.selfloop_edges())
        return network

