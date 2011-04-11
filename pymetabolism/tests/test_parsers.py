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


import nose.tools as nt

from pymetabolism.parsers import SBMLParser


class TestSBMLParser(object):

    def setUp(self):
        self.parser = SBMLParser()
        self.parser.parse("data/Ec_core_flux1.xml")

    def test_get_reactions(self):
        for rxn in self.parser.get_reactions():
            print rxn

    def test_get_compounds(self):
        for cmpd in self.parser.get_compounds():
            print cmpd

