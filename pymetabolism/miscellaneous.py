#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===================
Library Miscellanea
===================

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-07-01
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    miscellaneous.py
"""


import logging

from .singletonmixin import Singleton


class NullHandler(logging.Handler):
    """
    A stub logging handler that ignores all logging messages. This is the
    default handler for all library loggers.
    """

    def emit(self, record):
        pass


class OptionsManager(Singleton):
    """
    A centralised instance to handle some common options throughout this
    library.
    """

    def __init__(self, *args, **kw_args):
        Singleton.__init__(self, *args, **kw_args)
        self.compound_prefix = "M_"
        self.reaction_prefix = "R_"
        self.reversible_suffix = "_Rev"
        self.compartments = {"_c": "Cytosol", "_e": "Extra_organism",
                "_b": "Exchange", "_p": "Periplasm"}
        self.parser = "SBML"
        self.lp_solver = "gurobi"
        self.lower_bound = 0.0
        self.upper_bound = 1000.0
        self.numeric_threshold = 1E-09
        self.num_proc = 1

    def get_parser(self):
        """
        Returns
        -------
        A parser of the correct variant for the current parser attribute type.
        """
        from . import parsers
        if self.parser.lower() == "sbml":
            return parsers.SBMLParser()
        else:
            raise NotImplementedError("No other parsers are supported at"\
            " the moment")


