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


import nose.tools as nt

from ..miscellaneous import OptionsManager
from ..parsers import SBMLParser


class TestOptionsManager(object):

    def __init__(self):
        object.__init__(self)
        self.obj = OptionsManager.get_instance()

    def test_get_parser(self):
        parser = self.obj.get_parser()
        nt.assert_true(isinstance(parser, SBMLParser))
        self.obj.parser = "snafu"
        nt.assert_raises(NotImplementedError, self.obj.get_parser)

