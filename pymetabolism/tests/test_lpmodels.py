#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===============
LP Models Tests
===============

:Authors:
    Moritz Emanuel Beber
    Nikolaus Sonnenschein
:Date:
    2011-04-10
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    test_lpmodels.py
"""


import nose.tools as nt

from .. import lpmodels


class TestGurobiFacade(object):

    def __init__(self):
        self.options = OptionsManager()
        self.options.reversible_suffix = "r"
        self.parser = self.options.get_parser()
        self.system = self.parser.parse(os.path.join(os.path.dirname(__file__),
                "data", "Ec_core_flux1.xml"))
