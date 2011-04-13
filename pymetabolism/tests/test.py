#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys

from os import path


def run(verbosity=1, doctest=False):
    """Run pymetabolism tests.

    Parameters
    ----------
    verbosity: integer, optional
      Level of detail in test reports.  Higher numbers provide  more detail.

    doctest: bool, optional
      True to run doctests in code modules
    """
    try:
        import nose
    except ImportError:
        raise ImportError(
                "The nose package is needed to run the pymetabolism tests.")
    sys.stderr.write("Running pymetabolism tests:")
    install_dir = path.join(path.dirname(__file__), path.pardir)
    argv = [" ", "--verbosity=%d" % verbosity, "-w", install_dir, "-exe"]
    if doctest:
        argv.extend(["--with-doctest", "--doctest-extension=txt"])
    nose.run(argv=argv)

if __name__ == "__main__":
    sys.path.append(path.abspath("../.."))
    run()

