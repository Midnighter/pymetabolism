#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===========================
Flux Balance Analysis Model
===========================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
    Nikolaus Sonnenschein
:Date:
    2011-03-28
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    test_fba.py
"""


import os
import pymetabolism
import nose.tools as nt

from glob import glob


def setup(module):
    files = glob(os.path.join(os.path.dirname(__file__), "data", "*.xml"))
    module.systems = [pymetabolism.parse(filename) for filename in files]

def teardown(module):
    del module.systems

def test_fba():
    options = pymetabolism.OptionsManager.get_instance()
    for solver in ["gurobi"]:
        options.solver = solver
        yield TestFBAModel()

# need module level setup for OptionsManager, different models

class TestFBAModel(object):

    def __init__(self):
        self.options = pymetabolism.OptionsManager.get_instance()
        self.options.reversible_suffix = "r"
        self.parser = self.options.get_parser()
        self.system = self.parser.parse()
        (self.model, self.known_fluxes) = self.system.get_lp_model()

#    def test_knock_out_reaction(self):
#        """Tests if a reaction is correctly constrained to zero flux"""
#        glp = self.glp
#        glp.knock_out_reaction(['R("R_PGK_Rev")'])
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(), 0.5836778206)
#        glp.undo()
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(), 0.9259122)
#
#    def test_knock_out_metabolite(self):
#        """Tests if a metabolite is correctly constrained to zero"""
#        glp = self.glp
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(), 0.9259122)
#        glp.knock_out_metabolite(['Matpc', 'Madpc'])
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(), 0.)
#        glp.undo()
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(), 0.9259122)
#
#    def test_modify_reaction_bounds(self):
#        """Tests if reaction bounds are correctly modified."""
#        glp = self.glp
#        glp.modify_reaction_bounds({'R("R_PGK")':(0, 0), 'R("R_PGK_Rev")':(0, 0)})
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(), 0.5836778206)
#        glp.undo()
#        glp.undo()
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(), 0.9259122)
#
#    def test_modify_metabolite_bounds(self):
#        """Tests if metabolite bounds are correctly modified."""
#        glp=self.glp
#        glp.modify_metabolite_bounds({'Matpc':(0,0),'Madpc':(0,0)})
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(),0.)
#        glp.undo()
#        glp.undo()
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(),0.9259122)
#
#    def test_get_substrates_and_products(self):
#        """Tests if correct substrates and products are returned."""
#        (substrates, products) = self.glp.get_substrates_and_products('R("R_PGK")')
#        self.assertEqual(substrates, ('M3pgc', 'Matpc'))
#        self.assertEqual(products, ('M13dpgc', 'Madpc'))
#
#    def test_free_metabolites(self):
#        """Tests if metabolite constraints are correctly removed."""
#        glp=self.glp
#        glp.modify_metabolite_bounds({'Matpc':(0,0),'Madpc':(0,0)})
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(),0.)
#        glp.free_metabolites()
#        glp.simplex()
#        self.assertAlmostEqual(glp.getObjVal(),0.9259122)
#
#    def test_add_metabolite_drains(self):
#        """Tests if a drain for a certain metabolite is correctly added."""
#        glp=self.glp
#        glp.add_metabolite_drains('M3pgc')
#        temp=glp.get_column_names()
#        if 'M3pgs_Drain' in temp:
#            return True
#
#    def test_add_transporters(self):
#        """Tests if a transporter for a certain metabolite is correctly added."""
#        glp=self.glp
#        glp.add_trasporters('M3pgc')
#        temp=glp.get_column_names()
#        if 'M3pgs_Transp' in temp:
#            return True
#
#    def test_get_reaction_objective(self):
#        """Tests if the correct reaction objective is returned."""
#        obj = self.glp.get_reaction_objective()
#        self.assertEqual(obj, 'R("R_BiomassEcoli")')
#
#    def test_set_reaction_objective(self):
#        """Tests if the reaction objective is correctly set."""
#        obj = 'R("R_PGK")'
#        self.glp.set_objective(obj)
#        self.assertEqual(self.glp.get_reaction_objective(),obj)
#
