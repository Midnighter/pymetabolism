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


import nose.tools as nt
import pymetabolism.metabolism as pymet
import itertools


class TestBasicMetabolicComponent:

    def setUp(self):
        self.comp_1 = pymet.BasicMetabolicComponent("foo")
        self.comp_2 = pymet.BasicMetabolicComponent("bar")
        self.comp_3 = pymet.BasicMetabolicComponent()

    def test__str__(self):
        nt.assert_equal(str(self.comp_1), "foo")
        nt.assert_equal(str(self.comp_2), "bar")
        nt.assert_equal(str(self.comp_3), "BasicMetabolicComponent_%d" %
                self.comp_3._index)

    def test__repr__(self):
        nt.assert_equal(repr(self.comp_1),
                "<pymetabolism.metabolism.BasicMetabolicComponent, 1>")
        nt.assert_equal(repr(self.comp_2),
                "<pymetabolism.metabolism.BasicMetabolicComponent, 2>")
        nt.assert_equal(repr(self.comp_3),
                "<pymetabolism.metabolism.BasicMetabolicComponent, %d>" %
                self.comp_3._index)

    def test__lt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))

    def test__le__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))

    def test__eq__(self):
        nt.assert_equal(self.comp_1, pymet.BasicMetabolicComponent("foo"))
        nt.assert_equal(self.comp_2, pymet.BasicMetabolicComponent("bar"))
        nt.assert_equal(self.comp_3,
                pymet.BasicMetabolicComponent("BasicMetabolicComponent_%d" %
                self.comp_3._index))

    def test__ne__(self):
        for (x, y) in itertools.combinations([self.comp_1, self.comp_2,
                self.comp_3], 2):
            nt.assert_not_equal(x, y)

    def test__gt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))

    def test__ge__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))

    def test__cmp__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))


class TestBasicCompound:

    def setUp(self):
        self.comp_1 = pymet.BasicCompound("foo")
        self.comp_2 = pymet.BasicCompound("bar")
        self.comp_3 = pymet.BasicCompound()

    def test__str__(self):
        nt.assert_equal(str(self.comp_1), "foo")
        nt.assert_equal(str(self.comp_2), "bar")
        nt.assert_equal(str(self.comp_3), "BasicCompound_%d" %
                self.comp_3._index)

    def test__repr__(self):
        nt.assert_equal(repr(self.comp_1),
                "<pymetabolism.metabolism.BasicCompound, %d>" %
                self.comp_1._index)
        nt.assert_equal(repr(self.comp_2),
                "<pymetabolism.metabolism.BasicCompound, %d>" %
                self.comp_2._index)
        nt.assert_equal(repr(self.comp_3),
                "<pymetabolism.metabolism.BasicCompound, %d>" %
                self.comp_3._index)

    def test__lt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))

    def test__le__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))

    def test__eq__(self):
        nt.assert_equal(self.comp_1, pymet.BasicCompound("foo"))
        nt.assert_equal(self.comp_2, pymet.BasicCompound("bar"))
        nt.assert_equal(self.comp_3,
                pymet.BasicCompound("BasicCompound_%d" %
                self.comp_3._index))

    def test__ne__(self):
        for (x, y) in itertools.combinations([self.comp_1, self.comp_2,
                self.comp_3], 2):
            nt.assert_not_equal(x, y)

    def test__gt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))

    def test__ge__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))

    def test__cmp__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))


class TestBasicReaction:

    def setUp(self):
        self.comp_1 = pymet.BasicReaction("foo")
        self.comp_2 = pymet.BasicReaction("bar", True)
        self.comp_3 = pymet.BasicReaction()

    def test__str__(self):
        nt.assert_equal(str(self.comp_1), "foo")
        nt.assert_equal(str(self.comp_2), "bar")
        nt.assert_equal(str(self.comp_3), "BasicReaction_%d" %
                self.comp_3._index)

    def test__repr__(self):
        nt.assert_equal(repr(self.comp_1),
                "<pymetabolism.metabolism.BasicReaction, %d>" %
                self.comp_1._index)
        nt.assert_equal(repr(self.comp_2),
                "<pymetabolism.metabolism.BasicReaction, %d>" %
                self.comp_2._index)
        nt.assert_equal(repr(self.comp_3),
                "<pymetabolism.metabolism.BasicReaction, %d>" %
                self.comp_3._index)

    def test__lt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))

    def test__le__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))

    def test__eq__(self):
        nt.assert_equal(self.comp_1, pymet.BasicReaction("foo"))
        nt.assert_equal(self.comp_2, pymet.BasicReaction("bar"))
        nt.assert_equal(self.comp_3,
                pymet.BasicReaction("BasicReaction_%d" %
                self.comp_3._index))

    def test__ne__(self):
        for (x, y) in itertools.combinations([self.comp_1, self.comp_2,
                self.comp_3], 2):
            nt.assert_not_equal(x, y)

    def test__gt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))

    def test__ge__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))

    def test__cmp__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))

    def test_reversibility(self):
        nt.assert_false(self.comp_1.reversible)
        nt.assert_true(self.comp_2.reversible)
        nt.assert_equal(self.comp_2.reversible,
                pymet.BasicReaction("bar").reversible)


class TestSBMLCompartment:

    def setUp(self):
        self.comp_1 = pymet.SBMLCompartment("foo")
        self.comp_2 = pymet.SBMLCompartment("bar", True)

    def test__str__(self):
        nt.assert_equal(str(self.comp_1), "foo")
        nt.assert_equal(str(self.comp_2), "bar")

    def test__repr__(self):
        nt.assert_equal(repr(self.comp_1),
                "<pymetabolism.metabolism.SBMLCompartment, %d>" %
                self.comp_1._index)
        nt.assert_equal(repr(self.comp_2),
                "<pymetabolism.metabolism.SBMLCompartment, %d>" %
                self.comp_2._index)

    def test__lt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))

    def test__le__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))

    def test__eq__(self):
        nt.assert_equal(self.comp_1, pymet.SBMLCompartment("foo"))
        nt.assert_equal(self.comp_2, pymet.SBMLCompartment("bar"))

    def test__ne__(self):
        for (x, y) in itertools.combinations([self.comp_1, self.comp_2], 2):
            nt.assert_not_equal(x, y)

    def test__gt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))

    def test__ge__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))

    def test__cmp__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))

class TestSBMLCompound:

    def setUp(self):
        self.comp_1 = pymet.SBMLCompound("foo")
        self.comp_2 = pymet.SBMLCompound("bar")

    def test__str__(self):
        nt.assert_equal(str(self.comp_1), "foo")
        nt.assert_equal(str(self.comp_2), "bar")

    def test__repr__(self):
        nt.assert_equal(repr(self.comp_1),
                "<pymetabolism.metabolism.SBMLCompound, %d>" %
                self.comp_1._index)
        nt.assert_equal(repr(self.comp_2),
                "<pymetabolism.metabolism.SBMLCompound, %d>" %
                self.comp_2._index)

    def test__lt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))

    def test__le__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))

    def test__eq__(self):
        nt.assert_equal(self.comp_1, pymet.SBMLCompound("foo"))
        nt.assert_equal(self.comp_2, pymet.SBMLCompound("bar"))

    def test__ne__(self):
        for (x, y) in itertools.combinations([self.comp_1, self.comp_2], 2):
            nt.assert_not_equal(x, y)

    def test__gt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))

    def test__ge__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))

    def test__cmp__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))

    def test__contains__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__contains__, ("C",))


class TestSBMLCompartmentCompound:

    def setUp(self):
        self.comp_1 = pymet.SBMLCompartmentCompound(pymet.SBMLCompound("foo"),
                pymet.SBMLCompartment("bar"))
        self.comp_2 = pymet.SBMLCompartmentCompound(pymet.SBMLCompound("crow"),
                pymet.SBMLCompartment("bar"))

    def test__str__(self):
        nt.assert_equal(str(self.comp_1), "foo(bar)")
        nt.assert_equal(str(self.comp_2), "crow(bar)")

    def test__repr__(self):
        nt.assert_equal(repr(self.comp_1),
                "<pymetabolism.metabolism.SBMLCompartmentCompound, %d>" %
                self.comp_1._index)
        nt.assert_equal(repr(self.comp_2),
                "<pymetabolism.metabolism.SBMLCompartmentCompound, %d>" %
                self.comp_2._index)

    def test__lt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))

    def test__le__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))

    def test__eq__(self):
        nt.assert_equal(self.comp_1, pymet.SBMLCompartmentCompound(
                pymet.SBMLCompound("foo"), pymet.SBMLCompartment("bar")))
        nt.assert_equal(self.comp_2, pymet.SBMLCompartmentCompound(
                pymet.SBMLCompound("crow"), pymet.SBMLCompartment("bar")))

    def test__ne__(self):
        for (x, y) in itertools.combinations([self.comp_1, self.comp_2], 2):
            nt.assert_not_equal(x, y)

    def test__gt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))

    def test__ge__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))

    def test__cmp__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))

    def test__contains__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__contains__, ("C",))


class TestSBMLReaction:

    def setUp(self):
        self.comp_1 = pymet.SBMLReaction("foo", {pymet.SBMLCompound("A"): 2,
                pymet.SBMLCompound("B"): 1}, {pymet.SBMLCompound("C"): 1})
        self.comp_2 = pymet.SBMLReaction("bar", {pymet.SBMLCompound("D"): 2,
                pymet.SBMLCompound("E"): 1}, {pymet.SBMLCompound("F"): 1}, True)
        self.comp_3 = pymet.SBMLReaction("snafu",
                {pymet.SBMLCompartmentCompound(pymet.SBMLCompound("X"),
                pymet.SBMLCompartment("cyt")): 2,
                pymet.SBMLCompartmentCompound(pymet.SBMLCompound("Y"),
                pymet.SBMLCompartment("cyt")): 1},
                {pymet.SBMLCompartmentCompound(pymet.SBMLCompound("Z"),
                pymet.SBMLCompartment("cyt")): 1})

#    def test__str__(self):
#        nt.assert_equal(str(self.comp_1), "foo")
#        nt.assert_equal(str(self.comp_2), "bar")
#        nt.assert_equal(str(self.comp_3), "SBMLReaction_%d" %
#                self.comp_3._index)

    def test__repr__(self):
        nt.assert_equal(repr(self.comp_1),
                "<pymetabolism.metabolism.SBMLReaction, %d>" %
                self.comp_1._index)
        nt.assert_equal(repr(self.comp_2),
                "<pymetabolism.metabolism.SBMLReaction, %d>" %
                self.comp_2._index)
        nt.assert_equal(repr(self.comp_3),
                "<pymetabolism.metabolism.SBMLReaction, %d>" %
                self.comp_3._index)

    def test__lt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))

    def test__le__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))

    def test__eq__(self):
        nt.assert_equal(self.comp_1, pymet.SBMLReaction("foo",
                {pymet.SBMLCompound("A"): 2, pymet.SBMLCompound("B"): 1},
                {pymet.SBMLCompound("C"): 1}))
        nt.assert_equal(self.comp_2, pymet.SBMLReaction("bar",
                {pymet.SBMLCompound("D"): 2, pymet.SBMLCompound("E"): 1},
                {pymet.SBMLCompound("F"): 1}, True))
        nt.assert_equal(self.comp_3, pymet.SBMLReaction("snafu",
                {pymet.SBMLCompartmentCompound(pymet.SBMLCompound("X"),
                pymet.SBMLCompartment("cyt")): 2,
                pymet.SBMLCompartmentCompound(pymet.SBMLCompound("Y"),
                pymet.SBMLCompartment("cyt")): 1},
                {pymet.SBMLCompartmentCompound(pymet.SBMLCompound("Z"),
                pymet.SBMLCompartment("cyt")): 1}))

    def test__ne__(self):
        for (x, y) in itertools.combinations([self.comp_1, self.comp_2,
                self.comp_3], 2):
            nt.assert_not_equal(x, y)

    def test__gt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__gt__, (self.comp_2,))

    def test__ge__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__ge__, (self.comp_2,))

    def test__cmp__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__cmp__, (self.comp_2,))

    def test_reversibility(self):
        nt.assert_false(self.comp_1.reversible)
        nt.assert_true(self.comp_2.reversible)
        nt.assert_equal(self.comp_3.reversible, self.comp_1.reversible)

