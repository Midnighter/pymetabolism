#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=============
Parsers Tests
=============

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
    Nikolaus Sonnenschein
:Date:
    2011-04-10
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    test_parsers.py
"""


import os
import nose.tools as nt

from pymetabolism.parsers import SBMLParser


class TestSBMLParser(object):

    def setup(self):
        self.parser = SBMLParser()
        self.parser.parse(os.path.join(os.path.dirname(__file__), "data", "Ec_core_flux1.xml"))

    def test_compartments(self):
        for comp in self.parser.compartments:
            print comp

    def test_compounds(self):
        for cmpd in self.parser.compounds:
            print cmpd

    def test_reactions(self):
        for rxn in self.parser.reactions:
            print rxn

