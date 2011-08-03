#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=============
Parsers Tests
=============

:Authors:
    Moritz Emanuel Beber
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

from ..miscellaneous import OptionsManager
from .. import parsers
from ..metabolism import metabolism as pymet


class TestSBMLParser(object):

    def __init__(self):
        self.options = OptionsManager()
        self.options.reversible_suffix = "r"
        self.parser = parsers.SBMLParser()
        self.system = self.parser.parse(os.path.join(os.path.dirname(__file__),
                "data", "Ec_core_flux1.xml"))

    def test_compartments(self):
        nt.assert_equal(len(self.system.compartments), 3)
        for comp in self.system.compartments:
            nt.assert_true(isinstance(comp, pymet.SBMLCompartment))
            nt.assert_true(comp.name)

    def test_compounds(self):
        nt.assert_equal(len(self.system.compounds), 77)
        for cmpd in self.system.compounds:
            nt.assert_true(isinstance(cmpd, pymet.SBMLCompartmentCompound))
            nt.assert_true(cmpd.name)
            nt.assert_true(cmpd.identifier)
            nt.assert_true(len(cmpd.extended_name) > 0)

    def test_reactions(self):
        nt.assert_equal(len(self.system.reactions), 77)
        for rxn in self.system.reactions:
            nt.assert_true(isinstance(rxn, pymet.SBMLReaction))
            nt.assert_true(rxn.name)
            nt.assert_true(rxn.identifier)
            nt.assert_true(len(rxn.extended_name) > 0)

