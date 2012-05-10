#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=================================
Metabolic Network Representations
=================================

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-04-13
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    networks.py
"""


__all__ = ["CompoundCentricNetwork", "CompoundCentricMultiNetwork",
        "ReactionCentricNetwork", "ReactionCentricMultiNetwork",
        "MetabolicNetwork"]

import logging
import re
import itertools
import networkx as nx

from ..metabolism import metabolism as pymet
from ..errors import PyMetabolismError
from .. import miscellaneous as misc


logger = logging.getLogger(__name__)
logger.addHandler(misc.NullHandler())

options = misc.OptionsManager.get_instance()


class CompoundCentricNetwork(nx.DiGraph):
    """
    """

    def __init__(self, name="", *args, **kw_args):
        """
        """
        nx.DiGraph.__init__(self, name=name, *args, **kw_args)

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=True)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add compound nodes
        for (i, cmpd) in enumerate(self.nodes_iter()):
            indeces[cmpd] = i
            net.add_node(i, label=cmpd.name, shape="circle", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


class CompoundCentricMultiNetwork(nx.MultiDiGraph):
    """
    """

    def __init__(self, name="", *args, **kw_args):
        """
        """
        nx.MultiDiGraph.__init__(self, name=name, *args, **kw_args)

    def to_directed(self):
        """
        Return a copy with no multiple edges and no attributes.
        """
        copy = CompoundCentricNetwork(name="directed " + self.name)
        copy.add_nodes_from(self.nodes_iter())
        copy.add_edges_from(self.edges_iter())
        return copy

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=False)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add compound nodes
        for (i, cmpd) in enumerate(self.nodes_iter()):
            indeces[cmpd] = i
            net.add_node(i, label=cmpd.name, shape="circle", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


class ReactionCentricNetwork(nx.DiGraph):
    """
    """

    def __init__(self, name="", *args, **kw_args):
        """
        """
        nx.DiGraph.__init__(self, name=name, *args, **kw_args)

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        """
        """
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=False)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add reaction nodes
        for (i, rxn) in enumerate(self.nodes_iter()):
            indeces[rxn] = i
            net.add_node(i, label=rxn.name, shape="box", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


class ReactionCentricMultiNetwork(nx.MultiDiGraph):
    """
    """

    def __init__(self, name="", *args, **kw_args):
        """
        """
        nx.MultiDiGraph.__init__(self, name=name, *args, **kw_args)

    def to_directed(self):
        """
        Return a copy with no multiple edges and no attributes.
        """
        copy = CompoundCentricNetwork(name="directed_" + self.name)
        copy.add_nodes_from(self.nodes_iter())
        copy.add_edges_from(self.edges_iter())
        return copy

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=False)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add reaction nodes
        for (i, rxn) in enumerate(self.nodes_iter()):
            indeces[rxn] = i
            net.add_node(i, label=rxn.name, shape="box", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


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
            raise PyMetabolismError("unidentified metabolic component")
        if isinstance(v, pymet.BasicReaction):
            self.reactions.add(v)
        elif isinstance(v, pymet.BasicCompound):
            self.compounds.add(v)
        else:
            raise PyMetabolismError("unidentified metabolic component")
        nx.DiGraph.add_edge(self, u, v, **kw_args)

    def remove_edge(self, u, v):
        """
        """
        if isinstance(u, pymet.BasicReaction):
            self.reactions.remove(u)
        elif isinstance(u, pymet.BasicCompound):
            self.compounds.remove(u)
        else:
            raise PyMetabolismError("unidentified metabolic component")
        if isinstance(v, pymet.BasicReaction):
            self.reactions.remove(v)
        elif isinstance(v, pymet.BasicCompound):
            self.compounds.remove(v)
        else:
            raise PyMetabolismError("unidentified metabolic component")
        nx.DiGraph.remove_edge(self, u, v)

    def add_node(self, n, **kw_args):
        """
        """
        if isinstance(n, pymet.BasicReaction):
            self.reactions.add(n, **kw_args)
        elif isinstance(n, pymet.BasicCompound):
            self.compounds.add(n, **kw_args)
        else:
            raise PyMetabolismError("unidentified metabolic component")
        nx.DiGraph.add_node(self, n)

    def remove_node(self, n):
        """
        """
        if isinstance(n, pymet.BasicReaction):
            self.reactions.remove(n)
        elif isinstance(n, pymet.BasicCompound):
            self.compounds.remove(n)
        else:
            raise PyMetabolismError("unidentified metabolic component")
        nx.DiGraph.remove_node(self, n)

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

    def read_kegg(self, organism, wsdl="http://soap.genome.jp/KEGG.wsdl",
            browse=10):
        """
        Requires SOAPpy and an active internet connection.
        """
        soap_py = misc.load_module("SOAPpy",
                url="http://pywebsvcs.sourceforge.net/")
        # establish connection to DBGET server
        serv = soap_py.WSDL.Proxy(wsdl)
        # find appropriate organism
        # require user choice here, too many options
        choices = serv.bfind("genome " + organism)
        choices = [choice for choice in choices.split("\n") if choice]
        start = 0
        end = start + browse
        length = len(choices)
        searching = True
        while start < length and searching:
            msg = [""]
            msg.append("Showing organisms %d-%d of %d, please choose an index:"\
                    % (start, end - 1, length))
            msg.append("")
            for i in range(start, end):
                msg.append("[%d] %s" % (i, choices[i]))
            if end < length:
                msg.append("")
                msg.append("Type any non-integer to show the next %d organisms."\
                        % min(length - end, browse))
            msg.append("")
            try:
                selection = int(raw_input("\n".join(msg)))
            except ValueError:
                start = end
                end = min(length, end + browse)
            else:
                while True:
                    try:
                        choice = choices[selection]
                    except IndexError:
                        try:
                            selection = int(raw_input("Chosen index is outside"\
                                    " the allowed range, try again:\n"))
                        except ValueError:
                            pass
                    else:
                        searching = False
                        break
        logger.info("Please be patient, this will take a few minutes.")
        pattern = re.compile(r"genome:T\d+ (\w+),")
        mobj = pattern.match(choice)
        organism = mobj.group(1)
        pathways = serv.list_pathways(organism)
        logger.info("KEGG contains %d pathways for this organism.",
                len(pathways))
        pattern = re.compile(r"\w+")
        reactions = set()
        for path in pathways:
            reactions.update(serv.get_reactions_by_pathway(path.entry_id))
        logger.info("The pathways contain %d unique reactions", len(reactions))
        for rxn in reactions:
            logger.debug(rxn)
            info = serv.bget(rxn)
            info = info.split("\n")
            begin = -1
            stop = 0
            for (i, line) in enumerate(info):
                if line.startswith("RPAIR"):
                    begin = i
                    continue
                if begin > -1 and pattern.match(line):
                    stop = i
                    break
            if begin < 0:
                logger.warn("No reaction pair information for '%s'.", rxn)
                continue
            pairs = [info[begin].split()[1:]]
            for i in range(begin + 1, stop):
                pairs.append(info[i].split())
            logger.debug(str(pairs))
            reac = pymet.BasicReaction(rxn[3:], reversible=True)
            for line in pairs:
                (u, v) = line[1].split("_")
                self.add_edge(pymet.BasicCompound(u), reac, rpair=line[2])
                self.add_edge(reac, pymet.BasicCompound(v), rpair=line[2])

    def read_edgelist(self, path, delimiter=None, comments="#"):
        """
        """

        def build_node(name):
            if name.startswith(options.compound_prefix):
                compound = pymet.BasicCompound(name[len(options.compound_prefix):])
                return compound
            elif name.startswith(options.reaction_prefix):
                if name.endswith(options.reversible_suffix):
                    reaction = pymet.BasicReaction(name[len(options.reaction_prefix):
                            -len(options.reversible_suffix)], reversible=True)
                    reaction.reversible = True
                else:
                    reaction = pymet.BasicReaction(name[len(options.reaction_prefix):])
                return reaction
            else:
                raise TypeError("unrecognised metabolic object")

        options = misc.OptionsManager.get_instance()
        with open(path, "r") as file_handle:
            lines = [line.strip() for line in file_handle]
        for line in lines:
            if line.startswith(comments) or line == "":
                continue
            tmp = line.split(delimiter)
            u = build_node(tmp[0])
            if isinstance(u, pymet.BasicReaction) and\
                    tmp[0].endswith(options.reversible_suffix):
                continue
            v = build_node(tmp[1])
            if isinstance(v, pymet.BasicReaction) and\
                    tmp[1].endswith(options.reversible_suffix):
                continue
            self.add_edge(u, v)

    def write_edgelist(self, path, distinct=True, delimiter="\t", comments="#"):
        """
        """
        options = misc.OptionsManager.get_instance()
        lines = list()
        for rxn in self.reactions:
            rxn_name = options.reaction_prefix + rxn.name
            if rxn.reversible:
                if distinct:
                    rev_name = "%s%s%s" % (options.reaction_prefix, rxn.name,
                            options.reversible_suffix)
                else:
                    rev_name = rxn_name
                for cmpd in self.successors_iter(rxn):
                    lines.append("%s%s%s\n" % (rxn_name, delimiter, options.compound_prefix
                            + cmpd.name))
                    lines.append("%s%s%s\n" % (options.compound_prefix + cmpd.name,
                            delimiter, rev_name))
                for cmpd in self.predecessors_iter(rxn):
                    lines.append("%s%s%s\n" % (options.compound_prefix + cmpd.name,
                            delimiter, rxn_name))
                    lines.append("%s%s%s\n" % (rev_name, delimiter, options.compound_prefix
                            + cmpd.name))
            else:
                for cmpd in self.successors_iter(rxn):
                    lines.append("%s%s%s\n" % (rxn_name, delimiter, options.compound_prefix
                            + cmpd.name))
                for cmpd in self.predecessors_iter(rxn):
                    lines.append("%s%s%s\n" % (options.compound_prefix + cmpd.name,
                            delimiter, rxn_name))
        with open(path, "w") as file_handle:
            file_handle.writelines(lines)

    def to_compound_centric(self):
        """
        """

        def add_bi(u, v):
            network.add_edge(u, v)
            network.add_edge(v, u)

        network = CompoundCentricMultiNetwork("compound_centric_" + self.name)
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
        network.remove_edges_from(network.selfloop_edges())
        return network

    def to_reaction_centric(self):
        """
        """
        network = ReactionCentricMultiNetwork("reaction_centric_" + self.name)
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

    def draw(self, filename, output_format="pdf", layout_program="fdp",
                layout_args="", distinct=False):
        import pygraphviz as pgv
        options = misc.OptionsManager.get_instance()
        net = pgv.AGraph(directed=True, name=filename, strict=True)
        node_attr= dict()
        link_attr= dict()
        # add compound nodes
        indeces = dict(itertools.izip(self.compounds, itertools.count()))
        for (cmpd, i) in indeces.iteritems():
            net.add_node(i, label=cmpd.name, shape="circle", **node_attr)
        # add reactions
        indeces.update(itertools.izip(self.reactions,
                itertools.count(len(self.compounds))))
        i = len(self.compounds) + len(self.reactions)
        for rxn in self.reactions:
            net.add_node(indeces[rxn], label=rxn.name, shape="box", **node_attr)
            # add forward reaction links
            for cmpd in self.predecessors(rxn):
                net.add_edge(indeces[cmpd], indeces[rxn], **link_attr)
            for cmpd in self.successors(rxn):
                net.add_edge(indeces[rxn], indeces[cmpd], **link_attr)
            if rxn.reversible:
                if distinct:
                    rev = pymet.BasicReaction(rxn.name + options.reversible_suffix)
                    indeces[rev] = i
                    net.add_node(i, label=rev.name, shape="box", **node_attr)
                    # add backward reaction links
                    for cmpd in self.predecessors(rxn):
                        net.add_edge(indeces[rev], indeces[cmpd], **link_attr)
                    for cmpd in self.successors(rxn):
                        net.add_edge(indeces[cmpd], indeces[rev], **link_attr)
                    i += 1
                else:
                    # add backward reaction links
                    for cmpd in self.predecessors(rxn):
                        net.add_edge(indeces[rxn], indeces[cmpd],
                                style="dotted", **link_attr)
                    for cmpd in self.successors(rxn):
                        net.add_edge(indeces[cmpd], indeces[rxn],
                                style="dotted", **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)

    def to_system(self):
        system = pymet.MetabolicSystem(name=self.name)
        for rxn in self.reactions:
            subs = dict((pymet.SBMLCompound(str(cmpd)),
                    self[cmpd][rxn]["coefficient"]) for cmpd in self.pred[rxn])
            prods = dict((pymet.SBMLCompound(str(cmpd)),
                    self[rxn][cmpd]["coefficient"]) for cmpd in self.succ[rxn])
            if rxn.reversible:
                system.add(pymet.SBMLReaction(str(rxn), subs, prods,
                        rxn.reversible, lower_bound=-options.upper_bound,
                        upper_bound=options.upper_bound))
            else:
                system.add(pymet.SBMLReaction(str(rxn), subs, prods,
                        rxn.reversible, lower_bound=options.lower_bound,
                        upper_bound=options.upper_bound))
        return system

