#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===================
Miscellaneous Tests
===================

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-08-02
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    test_miscellaneous.py
"""


import os
import nose.tools as nt

from ..miscellaneous import OptionsManager
from ..lpmodels import GurobiFacade
from ..parsers import SBMLParser


class TestOptionsManager(object):

    def __init__(self):
        object.__init__(self)
        self.obj = OptionsManager()

    def test_get_parser(self):
        parser = self.obj.get_parser()
        nt.assert_true(isinstance(parser, SBMLParser))
        self.obj.parser = "snafu"
        nt.assert_raises(NotImplementedError, self.obj.get_parser)

    def test_get_lp_model(self):
        model = self.obj.get_lp_model()
        nt.assert_true(isinstance(model, GurobiFacade))
        model = self.obj.get_lp_model("foomanchu")
        nt.assert_equal(model.name, "foomanchu")
        self.obj.lp_solver = "snafu"
        nt.assert_raises(NotImplementedError, self.obj.get_lp_model)

