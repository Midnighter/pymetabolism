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


class TestSimpleCompound:

    def setUp(self):
        self.comp_1 = pymet.SimpleCompound("foo")
        self.comp_2 = pymet.SimpleCompound("bar")
        self.comp_3 = pymet.SimpleCompound()

    def test__str__(self):
        nt.assert_equal(str(self.comp_1), "foo")
        nt.assert_equal(str(self.comp_2), "bar")
        nt.assert_equal(str(self.comp_3), "SimpleCompound_%d" %
                self.comp_3._index)

    def test__repr__(self):
        nt.assert_equal(repr(self.comp_1),
                "<pymetabolism.metabolism.SimpleCompound, %d>" %
                self.comp_1._index)
        nt.assert_equal(repr(self.comp_2),
                "<pymetabolism.metabolism.SimpleCompound, %d>" %
                self.comp_2._index)
        nt.assert_equal(repr(self.comp_3),
                "<pymetabolism.metabolism.SimpleCompound, %d>" %
                self.comp_3._index)

    def test__lt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))

    def test__le__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))

    def test__eq__(self):
        nt.assert_equal(self.comp_1, pymet.SimpleCompound("foo"))
        nt.assert_equal(self.comp_2, pymet.SimpleCompound("bar"))
        nt.assert_equal(self.comp_3,
                pymet.SimpleCompound("SimpleCompound_%d" %
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


class TestSimpleReaction:

    def setUp(self):
        self.comp_1 = pymet.SimpleReaction("foo")
        self.comp_2 = pymet.SimpleReaction("bar", True)
        self.comp_3 = pymet.SimpleReaction()

    def test__str__(self):
        nt.assert_equal(str(self.comp_1), "foo")
        nt.assert_equal(str(self.comp_2), "bar")
        nt.assert_equal(str(self.comp_3), "SimpleReaction_%d" %
                self.comp_3._index)

    def test__repr__(self):
        nt.assert_equal(repr(self.comp_1),
                "<pymetabolism.metabolism.SimpleReaction, %d>" %
                self.comp_1._index)
        nt.assert_equal(repr(self.comp_2),
                "<pymetabolism.metabolism.SimpleReaction, %d>" %
                self.comp_2._index)
        nt.assert_equal(repr(self.comp_3),
                "<pymetabolism.metabolism.SimpleReaction, %d>" %
                self.comp_3._index)

    def test__lt__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))

    def test__le__(self):
        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))

    def test__eq__(self):
        nt.assert_equal(self.comp_1, pymet.SimpleReaction("foo"))
        nt.assert_equal(self.comp_2, pymet.SimpleReaction("bar"))
        nt.assert_equal(self.comp_3,
                pymet.SimpleReaction("SimpleReaction_%d" %
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
                pymet.SimpleReaction("bar").reversible)


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
#                "<pymetabolism.metabolism.SBMLCompartmentCompound, %d>" %
#                self.comp_1._index)
#        nt.assert_equal(repr(self.comp_2),
#                "<pymetabolism.metabolism.SBMLCompartmentCompound, %d>" %
#                self.comp_2._index)
#
#    def test__lt__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__lt__, (self.comp_2,))
#
#    def test__le__(self):
#        nt.assert_raises(NotImplementedError, self.comp_1.__le__, (self.comp_2,))
#
#    def test__eq__(self):
#        nt.assert_equal(self.comp_1, pymet.SBMLCompartmentCompound(pymet.SBMLCompound("foo"),
#                pymet.SBMLCompartment("bar")))
#        nt.assert_equal(self.comp_2, pymet.SBMLCompartmentCompound(pymet.SBMLCompound("foo"),
#                pymet.SBMLCompartment("bar")))
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

    def test__str__(self):
        nt.assert_equal(str(self.comp_1), "foo")
        nt.assert_equal(str(self.comp_2), "bar")
        nt.assert_equal(str(self.comp_3), "SBMLReaction_%d" %
                self.comp_3._index)

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

#class TestSBMLCompound:
#
#    def setUp(self):
#        self.r1 = Reaction('ATPhyrdolysis', (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (-1,-1,1,1))
#        self.r2 = Reaction('GTPhyrdolysis', (Compound('gtp'), Compound('h2o')), (Compound('gdp'), Compound('pi')), (-1,-1,1,1))
#        self.r3 = Reaction('CTPhyrdolysis', (Compound('ctp'), Compound('h2o')), (Compound('cdp'), Compound('pi')), (-1,-1,1,1))
#        self.metabolic_system = Metabolism((self.r1, self.r2, self.r3))
#        self.smaller_metabolic_system = Metabolism((self.r1, self.r2))
#    
#    # def test__str__(self):
#    #     """Tests if the __str__ methods provides the correct statistics"""
#    #     self.assertTrue(type(self.metabolic_system.__str__()) == str)
#        
#    def test__getitem__(self):
#        """Test if metabolic_system[y] provides the correct x"""
#        self.assertEqual(self.metabolic_system['ATPhyrdolysis'], self.r1)
#        self.assertEqual(self.metabolic_system[0], self.r1)
#        self.assertEqual(self.metabolic_system['GTPhyrdolysis'], self.r2)
#        self.assertEqual(self.metabolic_system[1], self.r2)
#
#    def test__contains__(self):
#        """Tests if __contains__ works correctly"""
#        self.assertTrue(self.r1 in self.metabolic_system)
#        self.assertTrue(self.r2 in self.metabolic_system)
#        self.assertTrue(self.r3 in self.metabolic_system)
#        self.assertFalse(self.r3 in self.smaller_metabolic_system)
#        
#        self.assertTrue('ATPhyrdolysis' in self.metabolic_system)
#        self.assertTrue('GTPhyrdolysis' in self.metabolic_system)
#        self.assertTrue('CTPhyrdolysis' in self.metabolic_system)
#        self.assertFalse('CTPhyrdolysis' in self.smaller_metabolic_system)
#        
#    
#    
#    # def test__add__(self):
#    #     """Tests if it is possible to merge metabolic systems"""
#    #     pass
#    # 
#    # def test__del__(self):
#    #     """docstring for test__del__"""
#    #     pass
#    # 
#    # def test__cmp__(self):
#    #     """docstring for test__"""
#    #     pass
#
#
#
#class Item4ReactionTests(unittest.TestCase):
#    def setUp(self):
#        self.reaction = Reaction('ATPhyrdolysis', (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (1,1,1,1))
#    
#    # def test__init__(self):
#    #     """Tests if wrong initializations are caught by the right Exceptions"""
#    #     self.assertRaises(TypeError, self.reaction.__init__, ('asdf', (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (1,1,1,1)))
#    #     self.assertRaises(TypeError, self.reaction.__init__, ([1,2], (Compound('atp'), Compound('h2o')), (Compound('adp'), Compound('pi')), (1,1,1,1)))
#    
#    def test__str__(self):
#        """Tests if the __str__ method works correctly"""
#        self.assertEqual(self.reaction.__str__(), '1 atp + 1 h2o -> 1 adp + 1 pi')
#
#    def test_contains_string_input(self):
#        self.assertTrue('atp' in self.reaction)
#        
#    def test_contains_string_input(self):
#        self.assertTrue(Compound('atp') in self.reaction)
#
#
#class Item2CompartmentTests(unittest.TestCase):
#    def setUp(self):
#        args = {"suffix": "_c", "spatial_dimensions": 3, "size": 1., "units": "ml"}
#        self.comp = Compartment("Cytosol", True, args)
#
#    def test_options(self):
#        def utility():
#            self.comp.options = None
#        self.assertRaises(AttributeError, utility())
#
#
#class Item3CompartCompoundTests(unittest.TestCase):
#    pass
#
#
#class Item1CompoundTests(unittest.TestCase):
#    pass
#
#
#if __name__ == '__main__':
#    tests = list()
#    tests.append(unittest.TestLoader().loadTestsFromTestCase(Item1CompoundTests))
#    tests.append(unittest.TestLoader().loadTestsFromTestCase(Item2CompartmentTests))
#    tests.append(unittest.TestLoader().loadTestsFromTestCase(Item3CompartCompoundTests))
#    tests.append(unittest.TestLoader().loadTestsFromTestCase(Item4ReactionTests))
#    tests.append(unittest.TestLoader().loadTestsFromTestCase(Item5MetabolismTests))
#    suite = unittest.TestSuite(tests)
#    unittest.TextTestRunner(verbosity=4).run(suite)
