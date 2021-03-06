#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
==========================
Metabolic Components Tests
==========================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
    Nikolaus Sonnenschein
:Date:
    2011-04-08
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    test_metabolism.py
"""


import inspect
import itertools
import nose.tools as nt

from .. import metabolism as pymet


def test_components():
    for (name, cls) in inspect.getmembers(pymet, predicate=inspect.isclass):
        if name.startswith("Basic"):
            yield (check_basic_cls, cls)
        elif name.startswith("SBML"):
            yield (check_sbml_cls, cls)
        elif name == "MetabolicSystem":
            yield (check_system, cls)

def check_basic_cls(cls):
    names = ["foo", "bar", None]
    instances = [cls(name) for name in names]
    # test __str__
    check__str__(instances, names)
    # test __repr__
    check__repr__(instances)

def check__str__(instances, names):
    for (instance, name) in itertools.izip(instances, names):
        if name == None:
            nt.assert_equal(str(instance),
                    "%s_%d" % (instance.__class__.__name__, instance._index))
        else:
            nt.assert_equal(str(instance), name)

def check__repr__(instances):
    for instance in instances:
        cls = instance.__class__
        nt.assert_equal(repr(instance), "<%s.%s, %d>" % (cls.__module__,
            cls.__name__, id(instance)))

def check_sbml_cls(cls):
    if cls.__name__ == "SBMLCompartment":
        pass
    elif cls.__name__ == "SBMLCompound":
        pass
    elif cls.__name__ == "SBMLCompartmentCompound":
        pass
    elif cls.__name__ == "SBMLReaction":
        pass

def check_system(cls):
    pass


#class TestSBMLCompartment:
#
#    def setUp(self):
#        self.comp_1 = pymet.SBMLCompartment("foo")
#        self.comp_2 = pymet.SBMLCompartment("bar", True)
#
#    def test__str__(self):
#        nt.assert_equal(str(self.comp_1), "foo")
#        nt.assert_equal(str(self.comp_2), "bar")
#
#    def test__repr__(self):
#        nt.assert_equal(repr(self.comp_1),
#                "<pymetabolism.metabolism.metabolism.SBMLCompartment, %d>" %
#                self.comp_1._index)
#        nt.assert_equal(repr(self.comp_2),
#                "<pymetabolism.metabolism.metabolism.SBMLCompartment, %d>" %
#                self.comp_2._index)
#
#    def test__lt__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))
#
#    def test__le__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))
#
#    def test__eq__(self):
#        nt.assert_equal(self.comp_1, pymet.SBMLCompartment("foo"))
#        nt.assert_equal(self.comp_2, pymet.SBMLCompartment("bar"))
#
#    def test__ne__(self):
#        for (x, y) in itertools.combinations([self.comp_1, self.comp_2], 2):
#            nt.assert_not_equal(x, y)
#
#    def test__gt__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))
#
#    def test__ge__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))
#
#    def test__cmp__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))
#
#class TestSBMLCompound:
#
#    def setUp(self):
#        self.comp_1 = pymet.SBMLCompound("foo")
#        self.comp_2 = pymet.SBMLCompound("bar")
#
#    def test__str__(self):
#        nt.assert_equal(str(self.comp_1), "foo")
#        nt.assert_equal(str(self.comp_2), "bar")
#
#    def test__repr__(self):
#        nt.assert_equal(repr(self.comp_1),
#                "<pymetabolism.metabolism.metabolism.SBMLCompound, %d>" %
#                self.comp_1._index)
#        nt.assert_equal(repr(self.comp_2),
#                "<pymetabolism.metabolism.metabolism.SBMLCompound, %d>" %
#                self.comp_2._index)
#
#    def test__lt__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))
#
#    def test__le__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))
#
#    def test__eq__(self):
#        nt.assert_equal(self.comp_1, pymet.SBMLCompound("foo"))
#        nt.assert_equal(self.comp_2, pymet.SBMLCompound("bar"))
#
#    def test__ne__(self):
#        for (x, y) in itertools.combinations([self.comp_1, self.comp_2], 2):
#            nt.assert_not_equal(x, y)
#
#    def test__gt__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))
#
#    def test__ge__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))
#
#    def test__cmp__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))
#
#    def test__contains__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__contains__, ("C",))
#
#
#class TestSBMLCompartmentCompound:
#
#    def setUp(self):
#        self.comp_1 = pymet.SBMLCompartmentCompound(pymet.SBMLCompound("foo"),
#                pymet.SBMLCompartment("bar"))
#        self.comp_2 = pymet.SBMLCompartmentCompound(pymet.SBMLCompound("crow"),
#                pymet.SBMLCompartment("bar"))
#
#    def test__str__(self):
#        nt.assert_equal(str(self.comp_1), "foo(bar)")
#        nt.assert_equal(str(self.comp_2), "crow(bar)")
#
#    def test__repr__(self):
#        nt.assert_equal(repr(self.comp_1),
#                "<pymetabolism.metabolism.metabolism.SBMLCompartmentCompound, %d>" %
#                self.comp_1._index)
#        nt.assert_equal(repr(self.comp_2),
#                "<pymetabolism.metabolism.metabolism.SBMLCompartmentCompound, %d>" %
#                self.comp_2._index)
#
#    def test__lt__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))
#
#    def test__le__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))
#
#    def test__eq__(self):
#        nt.assert_equal(self.comp_1, pymet.SBMLCompartmentCompound(
#                pymet.SBMLCompound("foo"), pymet.SBMLCompartment("bar")))
#        nt.assert_equal(self.comp_2, pymet.SBMLCompartmentCompound(
#                pymet.SBMLCompound("crow"), pymet.SBMLCompartment("bar")))
#
#    def test__ne__(self):
#        for (x, y) in itertools.combinations([self.comp_1, self.comp_2], 2):
#            nt.assert_not_equal(x, y)
#
#    def test__gt__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))
#
#    def test__ge__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))
#
#    def test__cmp__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))
#
#    def test__contains__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__contains__, ("C",))
#
#
#class TestSBMLReaction:
#
#    def setUp(self):
#        self.comp_1 = pymet.SBMLReaction("foo", {pymet.SBMLCompound("A"): 2,
#                pymet.SBMLCompound("B"): 1}, {pymet.SBMLCompound("C"): 1})
#        self.comp_2 = pymet.SBMLReaction("bar", {pymet.SBMLCompound("D"): 2,
#                pymet.SBMLCompound("E"): 1}, {pymet.SBMLCompound("F"): 1}, True)
#        self.comp_3 = pymet.SBMLReaction("snafu",
#                {pymet.SBMLCompartmentCompound(pymet.SBMLCompound("X"),
#                pymet.SBMLCompartment("cyt")): 2,
#                pymet.SBMLCompartmentCompound(pymet.SBMLCompound("Y"),
#                pymet.SBMLCompartment("cyt")): 1},
#                {pymet.SBMLCompartmentCompound(pymet.SBMLCompound("Z"),
#                pymet.SBMLCompartment("cyt")): 1})
#
##    def test__str__(self):
##        nt.assert_equal(str(self.comp_1), "foo")
##        nt.assert_equal(str(self.comp_2), "bar")
##        nt.assert_equal(str(self.comp_3), "SBMLReaction_%d" %
##                self.comp_3._index)
#
#    def test__repr__(self):
#        nt.assert_equal(repr(self.comp_1),
#                "<pymetabolism.metabolism.metabolism.SBMLReaction, %d>" %
#                self.comp_1._index)
#        nt.assert_equal(repr(self.comp_2),
#                "<pymetabolism.metabolism.metabolism.SBMLReaction, %d>" %
#                self.comp_2._index)
#        nt.assert_equal(repr(self.comp_3),
#                "<pymetabolism.metabolism.metabolism.SBMLReaction, %d>" %
#                self.comp_3._index)
#
#    def test__lt__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))
#
#    def test__le__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))
#
#    def test__eq__(self):
#        nt.assert_equal(self.comp_1, pymet.SBMLReaction("foo",
#                {pymet.SBMLCompound("A"): 2, pymet.SBMLCompound("B"): 1},
#                {pymet.SBMLCompound("C"): 1}))
#        nt.assert_equal(self.comp_2, pymet.SBMLReaction("bar",
#                {pymet.SBMLCompound("D"): 2, pymet.SBMLCompound("E"): 1},
#                {pymet.SBMLCompound("F"): 1}, True))
#        nt.assert_equal(self.comp_3, pymet.SBMLReaction("snafu",
#                {pymet.SBMLCompartmentCompound(pymet.SBMLCompound("X"),
#                pymet.SBMLCompartment("cyt")): 2,
#                pymet.SBMLCompartmentCompound(pymet.SBMLCompound("Y"),
#                pymet.SBMLCompartment("cyt")): 1},
#                {pymet.SBMLCompartmentCompound(pymet.SBMLCompound("Z"),
#                pymet.SBMLCompartment("cyt")): 1}))
#
#    def test__ne__(self):
#        for (x, y) in itertools.combinations([self.comp_1, self.comp_2,
#                self.comp_3], 2):
#            nt.assert_not_equal(x, y)
#
#    def test__gt__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))
#
#    def test__ge__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))
#
#    def test__cmp__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))
#
#    def test_reversibility(self):
#        nt.assert_false(self.comp_1.reversible)
#        nt.assert_true(self.comp_2.reversible)
#        nt.assert_equal(self.comp_3.reversible, self.comp_1.reversible)
#
#
#class TestMetabolicSystem(object):
#
#    def __init__(self):
#        self.options = OptionsManager()
#        self.options.reversible_suffix = "r"
#        self.parser = self.options.get_parser()
#        self.system = self.parser.parse(os.path.join(os.path.dirname(__file__),
#                "..", "..", "tests", "data", "Ec_core_flux1.xml"))
#        self.empty = None
#
#    def setup(self):
#        self.empty = pymet.MetabolicSystem()
#
#    def test_add(self):
#        c = pymet.BasicCompound()
#        self.empty.add(c)
#        sc = pymet.SBMLCompound("foo")
#        self.empty.add(sc)
#        sc2 = pymet.SBMLCompound("man")
#        self.empty.add(sc2)
#        co = pymet.SBMLCompartment("choo")
#        self.empty.add(co)
#        r = pymet.BasicReaction()
#        self.empty.add(r)
#        sr = pymet.SBMLReaction("stuff", substrates={sc2: 2}, products={sc: 1})
#        self.empty.add(sr)
#        nt.assert_equals(len(self.empty.compounds), 3)
#        nt.assert_equals(len(self.empty.reactions), 2)
#        nt.assert_equals(len(self.empty.compartments), 1)
#        self.empty.add(sc)
#        nt.assert_equals(len(self.empty.compounds), 3)
#        nt.assert_raises(PyMetabolismError, self.empty.add, 4)
#
#    def test_update(self):
#        c = pymet.BasicCompound()
#        sc = pymet.SBMLCompound("foo")
#        sc2 = pymet.SBMLCompound("man")
#        self.empty.update([c], type(c))
#        self.empty.update([sc, sc2], type(sc))
#        co = pymet.SBMLCompartment("choo")
#        self.empty.update([co], type(co))
#        r = pymet.BasicReaction()
#        sr = pymet.SBMLReaction("stuff", substrates={sc2: 2}, products={sc: 1})
#        self.empty.update([r], type(r))
#        self.empty.update([sr], type(sr))
#        nt.assert_equals(len(self.empty.compounds), 3)
#        nt.assert_equals(len(self.empty.reactions), 2)
#        nt.assert_equals(len(self.empty.compartments), 1)
#        self.empty.update([c], type(c))
#        nt.assert_equals(len(self.empty.compounds), 3)
#        nt.assert_raises(PyMetabolismError, self.empty.update, [4], type(4))
#
#    def test_verify_consistency(self):
#        nt.assert_false(self.system.verify_consistency())
#        r = pymet.SBMLReaction("Biomass_Ecoli_core_N__w_GAM_", {}, {})
#        self.system.reactions.remove(r)
#        nt.assert_true(self.system.verify_consistency())
#        self.system = self.parser.parse(os.path.join(os.path.dirname(__file__),
#                "..", "..", "tests", "data", "Ec_core_flux1.xml"))
#
#    def test_detect_unconserved_metabolites(self):
#        pass
#
#    def test_generate_fba_model(self):
#        pass
#
#    def test_generate_network(self):
#        network = self.system.generate_network(stoichiometric_coefficients=True)
#        nt.assert_equals(network.order(), len(self.system.compounds) +
#                len(self.system.reactions))
#        nt.assert_equals(network.size(), sum(len(rxn) for rxn in\
#            self.system.reactions))
#
