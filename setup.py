#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
====================
PyMetabolism Package
====================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
    Nikolaus Sonnenschein
:Date:
    2011-04-12
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    setup.py
"""


from glob import glob
from distutils.core import setup


setup(
    name = "pymetabolism",
    version = "0.1",
    description = "handling and analysis of metabolic models",
    author = "Moritz Beber",
    author_email = "moritz (dot) beber (at) googlemail (dot) com",
    url = "http://github.com/Midnighter/pymetabolism",
    packages = ["pymetabolism",
            "pymetabolism.metabolism",
            "pymetabolism.network",
            "pymetabolism.tests",
            "pymetabolism.metabolism.tests",
            "pymetabolism.network.tests",],
    package_data = {"pymetabolism.tests": ["data/*xml", "data/*lp"]},
    )

