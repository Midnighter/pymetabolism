..  -*- coding: utf-8 -*-

.. currentmodule:: pymetabolism


Before starting any work with the pymetabolism package, you should be aware that
it uses logging facilities. Each submodule of pymetabolism logs certain
information that are muted by default. You can access that information by
creating a low-level logger and adding a handler to it. All information in
submodules will then be propagated to the handler of that base logger. For more
information see also the logging module documentation.

>>> import logging
>>> logger = logging.getLogger("pymetabolism")

You could also use the root logger logging.getLogger("").

>>> logger.setLevel(logging.INFO)
>>> logger.addHandler(logging.StreamHandler())

The next thing you need to be aware of is that pymetabolism uses a Singleton
instance to manage a few global options. There you can use some common suffixes
and prefixes used in parsing documents, printing information, and writing to
files.

With the OptionsManager you also determine what type of parser you want and the
type of solver you want to use for linear optimization (if any).

In many cases the starting point of using this package will be a metabolic model
in SBML format. That model contained in a file will have to be parsed. In the
following is a possible scenario of using the pymetabolism package.

>>> import pymetabolism
>>> options = pymetabolism.OptionsManager.get_instance()

For now the get_instance call is necessary due to the mechanics of the Singleton
class. The OptionsManager comes with a variety of reasonable default values in
place.

>>> print options.parser
'SBML'
>>> print options.lp_solver
'gurobi'

Let's say we have an SBML document containing a metabolic model interest. The
SBML format is the default so we just get a parser instance.

>>> parser = options.get_parser()

When parsing an SBML document a MetabolicSystem instance is returned.

>>> system = parser.parse("path/to/model.xml")

A MetabolicSystem is basically just a container for all the compartments,
compounds, and reactions found in the model. When you have a system you can
decide how to continue from there. A first step could be to check whether the
system adheres to mass conservation rules.

>>> system.verify_consistency()
True

From here on the system can be converted either into a linear programming model
suitable for flux balance analysis or into a bipartite network representing the
metabolism.

>>> model = system.generate_fba_model()
>>> model.fba()
2
>>> model.get_objective_value()
0.0

Please also consult the documentation for the FBAModel class.

The bipartite network is used using the networkx package and its DiGraph class.
So far few methods of the networkx.DiGraph class are overridden. If you do evil
things the bipartite nature may be destroyed and many networkx algorithms may
yield strange results due to the existing bipartite nature. So be careful.

>>> network = system.generate_network()

The MetabolicNetwork class can be further used to generate unipartite
projections and to draw it. Look at its documentation for further use cases.

