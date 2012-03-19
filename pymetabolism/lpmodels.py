#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
====================
LP Solver Interfaces
====================

:Authors:
    Moritz Emanuel Beber
    Nils Kölling
    Nikolaus Sonnenschein
:Date:
    2011-03-28
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    lpmodels.py
"""


import sys
import os
import itertools
import copy
import logging
import tempfile

from . import miscellaneous as misc
from .errors import PyMetabolismError


logger = logging.getLogger(__name__)
logger.addHandler(misc.NullHandler())


options = misc.OptionsManager.get_instance()


class MetaLPModelFacade(type):
    """
    Meta class that, according to the solver chosen in the options, populates
    the given class with methods that are solver specific.

    The general interface of the created class is supposed to look like the one
    described by `FBAModel`.
    """

    def __new__(mcls, name, bases, dct):
        """
        """
        if options.lp_solver.lower() == "gurobi":
            _grb_populate(dct)
        return super(MetaLPModelFacade, mcls).__new__(mcls, name, bases, dct)


################################################################################
# Gurobi Facade
################################################################################


def _grb_populate(attrs):
    name = "gurobipy"
    if not sys.modules.has_key(name):
        try:
            __import__(name)
        except ImportError:
            raise ImportError(u"Gurobi with python bindings is required for "\
                    "this functionality, please supply a different solver in "\
                    "the options or visit http://www.gurobi.com/ for "\
                    "detailed installation instructions.")

    grb = sys.modules[name]
    # suppress reports to stdout
    grb.setParam("OutputFlag", 0)
    # deal with Gurobi's annoying log file
    tmp_file = tempfile.mkstemp()[1] # absolute path component
    grb.setParam("LogFile", tmp_file)
    os.remove("gurobi.log")
    # set the number of processes

    # set class attributes
    attrs["_grb"] = grb
    for (key, value) in attrs.iteritems():
        if key.startswith("_"):
            continue
        try:
            attrs[key] = eval("_grb_" + key)
            attrs[key].__doc__ = value.__doc__
        except NameError:
            pass

    # gurobi helper functions
    attrs["_flux"] = _grb__flux
    attrs["_status"] = _grb__status
    attrs["_make_binary"] = _grb__make_binary
    attrs["_make_integer"] = _grb__make_integer

def _grb_initialise(self, name):
    self._model = self._grb.Model(name)
    self._rxn2var = dict()
    self._var2rxn = dict()
    self._rev2var = dict()
    self._var2rev = dict()
    self._cmpd2cnstrnt = dict()
    self._cnstrnt2cmpd = dict()
    self._sources = dict()
    self._drains = dict()

def _grb___copy__(self):
    # TODO
    cpy = FBAModel(self.name)
    cpy._model = self._model.copy()
    cpy._system2grb = dict()
    for col in cpy._model.getVars():
        cpy._system2grb[col.getAttr("VarName")] = col
    cpy._system2grb = dict()
    for row in cpy._model.getConstrs():
        cpy._system2grb[row.getAttr("ConstrName")] = row
    return cpy

def _grb___deepcopy__(self, memo=dict()):
    # TODO
    return self.__copy__()

def _grb__grb_copy(self):
    # TODO
    return self.__deepcopy__()

def _grb__make_binary(self, names):
    for col in names:
        var = self._variables[col]
        var.vType = "B"
    self._model.update()

def _grb__make_integer(self, names):
    for col in names:
        var = self._variables[col]
        var.vType = "I"
    self._model.update()

def _grb_add_reaction(self, reaction, coefficients, lb=None, ub=None):
    if lb is None:
        lb = options.lower_bound
    if ub is None:
        ub = options.upper_bound
    if hasattr(reaction, "__iter__"):
    # we really add multiple reactions
        if not hasattr(lb, "__iter__"):
            lb_iter = itertools.repeat(lb)
        else:
            lb_iter = lb
        if not hasattr(ub, "__iter__"):
            ub_iter = itertools.repeat(ub)
        else:
            ub_iter = ub
        for (rxn, lb, ub) in itertools.izip(reaction, lb_iter, ub_iter):
            if self._rxn2var.has_key(rxn):
                continue
            if rxn.reversible:
                if lb < 0:
                    var_rev = self._model.addVar(0.0, abs(lb), name=str(rxn) +
                            options.reversible_suffix)
                    var = self._model.addVar(0.0, ub, name=str(rxn))
                else:
                    var_rev = self._model.addVar(lb, ub, name=str(rxn) +
                            options.reversible_suffix)
                    var = self._model.addVar(lb, ub, name=str(rxn))
                self._rev2var[rxn] = var_rev
                self._var2rev[var_rev] = rxn
            else:
                var = self._model.addVar(lb, ub, name=str(rxn))
            self._rxn2var[rxn] = var
            self._var2rxn[var] = rxn
        self._model.update()
        for (rxn, coeff_iter) in itertools.izip(reaction, coefficients):
            # will throw KeyError if something went wrong above
            var = self._rxn2var[rxn]
            flag = False
            for (cmpd, factor) in coeff_iter:
                cnstrnt = self._cmpd2cnstrnt.get(cmpd)
                if cnstrnt:
                    self._model.chgCoeff(cnstrnt, var, factor)
                else:
                    flag = True
                    cnstrnt = self._model.addConstr(
                            self._grb.LinExpr(factor, var),
                            self._grb.GRB.EQUAL, 0.0, name=str(cmpd))
                    self._cmpd2cnstrnt[cmpd] = cnstrnt
            if flag:
                self._model.update()
            # add constraints with inverse factors for reverse reaction
            if rxn.reversible:
                var = self._rev2var[rxn]
                for (cmpd, factor) in coeff_iter:
                    cnstrnt = self._cmpd2cnstrnt[cmpd]
                    self._model.chgCoeff(cnstrnt, var, -factor)
    else:
        if self._rxn2var.has_key(reaction):
            return
        if reaction.reversible:
            if lb < 0:
                var_rev = self._model.addVar(0.0, abs(lb), name=str(reaction) +
                        options.reversible_suffix)
                var = self._model.addVar(0.0, ub, name=str(reaction))
            else:
                var_rev = self._model.addVar(lb, ub, name=str(reaction) +
                        options.reversible_suffix)
                var = self._model.addVar(lb, ub, name=str(reaction))
            self._rev2var[reaction] = var_rev
            self._var2rev[var_rev] = reaction
        else:
            var = self._model.addVar(lb, ub, name=str(reaction))
        self._rxn2var[reaction] = var
        self._var2rxn[var] = reaction
        self._model.update()
        flag = False
        for (cmpd, factor) in coefficients:
            cnstrnt = self._system2grb.get(cmpd)
            if cnstrnt:
                self._model.chgCoeff(cnstrnt, var, factor)
            else:
                flag = True
                cnstrnt = self._model.addConstr(
                        self._grb.LinExpr(factor, var),
                        self._grb.GRB.EQUAL, 0.0, name=str(cmpd))
                self._system2grb[cmpd] = cnstrnt
                self._grb2system[cnstrnt] = cmpd
        if flag:
            self._model.update()
        # add constraints with inverse factors for reverse reaction
        if reaction.reversible:
            var = self._rev2var[reaction]
            for (cmpd, factor) in coefficients:
                cnstrnt = self._cmpd2cnstrnt[cmpd]
                self._model.chgCoeff(cnstrnt, var, -factor)

def _grb_iter_reactions(self, compound=None, coefficients=False):
    if not compound:
        return self._rxn2var.iterkeys()
    result = list()
    if coefficients:
        lin_expr = self._model.getRow(self._cmpd2cnstrnt[compound])
        for i in range(lin_expr.size()):
            var = lin_expr.getVar(i)
            rxn = self._var2rxn[var] if self._var2rxn.has_key(var) else\
                    self._var2rev[var]
            factor = lin_expr.getCoeff(i)
            result.append((rxn, factor))
    else:
        for i in range(lin_expr.size()):
            var = lin_expr.getVar(i)
            rxn = self._var2rxn[var] if self._var2rxn.has_key(var) else\
                    self._var2rev[var]
            result.append(rxn)
    return result

def _grb_modify_reaction_bounds(self, reaction, lb=None, ub=None):
    # we allow for lazy updating of the model here (better not be a bug)
    if lb == None:
        lb = options.lower_bound
    if ub == None:
        ub = options.upper_bound
    if hasattr(reaction, "__iter__"):
        # we really modify multiple reactions
        if not hasattr(lb, "__iter__"):
            lb_iter = itertools.repeat(lb)
        if not hasattr(ub, "__iter__"):
            ub_iter = itertools.repeat(ub)
        for (rxn, lb, ub) in itertools.izip(reaction, lb_iter, ub_iter):
            var = self._rxn2var[rxn]
            if rxn.reversible:
                var_rev = self._rev2var[rxn]
                if lb < 0.0:
                    var_rev.lb = 0.0
                    var_rev.ub = abs(lb)
                    var.lb = 0.0
                    var.ub = ub
                else:
                    var_rev.lb = lb
                    var_rev.ub = ub
                    var.lb = lb
                    var.ub = ub
            else:
                var.lb = lb
                var.ub = ub
    else:
        var = self._system2grb[reaction]
        if reaction.reversible:
            var_rev = self._rev2var[reaction]
            if lb < 0.0:
                var_rev.lb = 0.0
                var_rev.ub = abs(lb)
                var.lb = 0.0
                var.ub = ub
            else:
                var_rev.lb = lb
                var_rev.ub = ub
                var.lb = lb
                var.ub = ub
        else:
            var.lb = lb
            var.ub = ub

def _grb_iter_reaction_bounds(self, reaction=None):
    # we rely on reversible reactions being treated in unison
    if reaction is None:
        reaction = self._rxn2var.iterkeys()
        return ((rxn, self._rxn2var[rxn].getAttr("LB"),
                self._rxn2var[rxn].getAttr("UB")) for rxn in reaction)
    elif hasattr(reaction, "__iter__"):
        # we really get multiple reactions
        return ((self._rxn2var[rxn].getAttr("LB"),
                self._rxn2var[rxn].getAttr("UB")) for rxn in reaction)
    else:
        var = self._rxn2var[reaction]
        return (var.lb, var.ub)

def _grb_modify_reaction_coefficients(self, reaction, coefficients):
    # we allow for lazy updating of the model here (better not be a bug)
    if hasattr(reaction, "__iter__"):
        for (rxn, coeff_iter) in itertools.izip(reaction, coefficients):
            var = self._rxn2var[rxn]
            for (cmpd, factor) in coeff_iter:
                self._model.chgCoeff(self._cmpd2cnstrnt[cmpd], var, factor)
            if rxn.reversible:
                var = self._rev2var[rxn]
                for (cmpd, factor) in coeff_iter:
                    self._model.chgCoeff(self._cmpd2cnstrnt[cmpd], var, -factor)
    else:
        var = self._rxn2var[reaction]
        for (cmpd, factor) in coefficients:
            self._model.chgCoeff(self._cmpd2cnstrnt[cmpd], var, factor)
        if reaction.reversible:
            var = self._rev2var[reaction]
            for (cmpd, factor) in coefficients:
                self._model.chgCoeff(self._cmpd2cnstrnt[cmpd], var, -factor)

def _grb_delete_reaction(self, reaction):
    if hasattr(reaction, "__iter__"):
        for rxn in reaction:
            var = self._rxn2var.pop(rxn)
            del self._var2rxn[var]
            self._model.remove(var)
            if rxn.reversible:
                var = self._rev2var.pop(rxn)
                del self._var2rev[var]
                self._model.remove(var)
    else:
        var = self._rxn2var.pop(reaction)
        del self._var2rxn[var]
        self._model.remove(var)
        if reaction.reversible:
            var = self._rev2var.pop(reaction)
            del self._var2rev[var]
            self._model.remove(var)
    self._model.update()

#def _grb_add_compound(self, compound, coefficients):
#    constraint = self._constraints.get(row)
#    if not constraint:
#        constraint = self._model.addConstr(0.0, self._gurobipy.GRB.EQUAL,
#                0.0, name=row)
#        self._constraints[row] = constraint
#        self._model.update()
#    for (column, factor) in coefficients:
#        var = self._variables.get(column)
#        if var:
#            self._model.chgCoeff(constraint, var, factor)
#        else:
#            raise PyMetabolismError("modifying coefficient of a"\
#                " non-existant column, please add the column first")
#    self._model.update()
#
#def _grb_add_rows(self, rows):
#    for row in rows:
#        if not row in self._constraints:
#            self._constraints[row] = self._model.addConstr(0.0,
#                    self._gurobipy.GRB.EQUAL, 0.0, name=row.name)
#    self._model.update()
#
#def _grb_delete_row(self, rows):
#    if hasattr(rows, "__iter__"):
#        for name in rows:
#            self._model.remove(self._constraints.pop(name))
#    else:
#        self._model.remove(self._variables.pop(rows))
#    self._model.update()
#
#def _grb_iter_compounds(self, column=None, coefficients=False):
#    if column:
#        grb_column = self._model.getCol(self._variables[column])
#        if coefficients:
#            return ((grb_column.getConstr(i).getAttr("ConstrName"),
#                grb_column.getCoeff(i)) for i in xrange(grb_column.size()))
#        else:
#            return (grb_column.getConstr(i).getAttr("ConstrName") for i in
#                xrange(grb_column.size()))
#    else:
#        return self._constraints.iterkeys()
#
#def _grb_modify_row_coefficients(self, row, coefficients):
#    constraint = self._constraints[row]
#    for (column, factor) in coefficients.iteritems():
#        self._model.chgCoeff(constraint, self._variables[column], factor)
#    self._model.update()
#
#def _grb_modify_row_bounds(self, rows):
#    raise NotImplementedError("Gurobi does not support bounds for rows")

def _grb_add_compound_source(self, compound, lb=None, ub=None):
    if lb is None:
        lb = options.lower_bound
    if ub is None:
        ub = options.upper_bound
    if hasattr(compound, "__iter__"):
    # we really add multiple compounds
        if hasattr(lb, "__iter__"):
            lb_iter = lb
        else:
            lb_iter = itertools.repeat(lb)
        if hasattr(ub, "__iter__"):
            ub_iter = ub
        else:
            ub_iter = itertools.repeat(ub)
        for (cmpd, lb, ub) in itertools.izip(compound, lb_iter, ub_iter):
            if self._sources.has_key[cmpd]:
                continue
            self._sources[cmpd] = self._model.addVar(lb, ub, name=str(cmpd) + "_Source")
        self._model.update()
        # we allow for lazy updating of the model here (better not be a bug)
        for cmpd in compound:
            var = self._sources[cmpd]
            cnstrnt = self._cmpd2cnstrnt[cmpd]
            self._model.chgCoeff(cnstrnt, var, 1.0)
    else:
        var = self._sources.get(compound)
        if not var:
            var = self._model.addVar(lb, ub, name=str(cmpd) + "_Source")
            self._sources[cmpd] = var
            self._model.update()
        cnstrnt = self._cmpd2cnstrnt[compound]
        self._model.chgCoeff(cnstrnt, var, 1.0)

def _grb_add_compound_drain(self, compound, lb=None, ub=None):
    if lb is None:
        lb = options.lower_bound
    if ub is None:
        ub = options.upper_bound
    if hasattr(compound, "__iter__"):
    # we really add multiple compounds
        if hasattr(lb, "__iter__"):
            lb_iter = lb
        else:
            lb_iter = itertools.repeat(lb)
        if hasattr(ub, "__iter__"):
            ub_iter = ub
        else:
            ub_iter = itertools.repeat(ub)
        for (cmpd, lb, ub) in itertools.izip(compound, lb_iter, ub_iter):
            if self._drains.has_key[cmpd]:
                continue
            self._drains[cmpd] = self._model.addVar(lb, ub, name=str(cmpd) +
                    "_Drain")
        self._model.update()
        # we allow for lazy updating of the model here (better not be a bug)
        for cmpd in compound:
            var = self._drains[cmpd]
            cnstrnt = self._cmpd2cnstrnt[cmpd]
            self._model.chgCoeff(cnstrnt, var, -1.0)
    else:
        var = self._drains.get(compound)
        if not var:
            var = self._model.addVar(lb, ub, name=str(cmpd) + "_Drain")
            self._drains[cmpd] = var
            self._model.update()
        cnstrnt = self._cmpd2cnstrnt[compound]
        self._model.chgCoeff(cnstrnt, var, -1.0)

def _grb_iter_objective_reaction(self, coefficients=False):
    lin_expr = self._model.getObjective()
    if coefficients:
        return ((self._var2rxn[lin_expr.getVar(i)], lin_expr.getCoeff(i))
                for i in range(lin_expr.size()))
    else:
        return (self._var2rxn[lin_expr.getVar(i)] for i in range(lin_expr.size()))

def _grb_set_objective(self, reaction):
    # we allow for lazy updating of the model here (better not be a bug)
    if hasattr(reaction, "__iter__"):
        for (rxn, factor) in reaction:
            self._rxn2var[rxn].setAttr("Obj", factor)
            if rxn.reversible:
                self._rev2var[rxn].setAttr("Obj", factor)
    else:
        self._rxn2var[reaction].setAttr("Obj", factor)
        if reaction.reversible:
            self._rev2var[reaction].setAttr("Obj", factor)

def _grb__status(self):
    """
    Determine the current status of the Gurobi model.
    """
    status = self._model.getAttr("Status")
    if status == self._gurobipy.GRB.LOADED:
        raise PyMetabolismError("optimize before retrieving the objective value", errorno=status)
    elif status == self._gurobipy.GRB.OPTIMAL:
        pass
    elif status == self._gurobipy.GRB.INFEASIBLE:
        raise PyMetabolismError("model is infeasible", errorno=status)
    elif status == self._gurobipy.GRB.INF_OR_UNBD:
        raise PyMetabolismError("model is infeasible or unbounded", errorno=status)
    elif status == self._gurobipy.GRB.UNBOUNDED:
        raise PyMetabolismError("model is unbounded", errorno=status)
    elif status == self._gurobipy.GRB.CUTOFF:
        raise PyMetabolismError("model solution is worse than provided cut-off", errorno=status)
    elif status == self._gurobipy.GRB.ITERATION_LIMIT:
        raise PyMetabolismError("iteration limit exceeded", errorno=status)
    elif status == self._gurobipy.GRB.NODE_LIMIT:
        raise PyMetabolismError("node limit exceeded", errorno=status)
    elif status == self._gurobipy.GRB.TIME_LIMIT:
        raise PyMetabolismError("time limit exceeded", errorno=status)
    elif status == self._gurobipy.GRB.SOLUTION_LIMIT:
        raise PyMetabolismError("solution limit reached", errorno=status)
    elif status == self._gurobipy.GRB.INTERRUPTED:
        raise PyMetabolismError("optimization process was interrupted", errorno=status)
    elif status == self._gurobipy.GRB.SUBOPTIMAL:
        raise PyMetabolismError("solution is suboptimal", errorno=status)
    elif status == self._gurobipy.GRB.NUMERIC:
        raise PyMetabolismError("optimization aborted due to numeric difficulties", errorno=status)

def _grb_get_objective_value(self):
    # _status should catch all problems (monitor this)
    self._status()
    return self._model.getAttr("ObjVal")

def _grb_fba(self, maximize=True):
    if maximize:
        self._model.setAttr("ModelSense", -1)
    else:
        self._model.setAttr("ModelSense", 1)
    self._model.optimize()

def _grb_parsimonious_fba(self, maximize=True):
    # implement minimization step and warm start
    if maximize:
        self._model.setAttr("ModelSense", -1)
    else:
        self._model.setAttr("ModelSense", 1)
    self._model.optimize()

def _grb__flux(self, reaction, threshold):
    flux = self._rxn2var[reaction].x
    if reaction.reversible:
        flux += self._rev2var[reaction].x
    return flux if flux > threshold else 0.0

def _grb_iter_flux(self, reaction=None, threshold=None):
    self._status()
    if threshold is None:
        threshold = options.numeric_threshold
    if reaction is None:
        reaction = self._rxn2var.iterkeys()
        return ((rxn, self._flux(rxn, threshold)) for rxn in reaction)
    elif hasattr(reaction, "__iter__"):
        return (self._flux(rxn, threshold) for rxn in reaction)
    else:
        return self._flux(reaction, threshold)

def _grb_set_medium(self, compound, lb=None, ub=None):
    # we allow for lazy updating of the model here (better not be a bug)
    if lb is None:
        lb = options.lower_bound
    if ub is None:
        ub = options.upper_bound
    # constrain all sources first
    for source in self._sources.itervalues():
        source.lb = 0.0
        source.ub = 0.0
    if hasattr(compound, "__iter__"):
    # we really add multiple compounds
        if hasattr(lb, "__iter__"):
            lb_iter = lb
        else:
            lb_iter = itertools.repeat(lb)
        if hasattr(ub, "__iter__"):
            ub_iter = ub
        else:
            ub_iter = itertools.repeat(ub)
        for (cmpd, lb, ub) in itertools.izip(compound, lb_iter, ub_iter):
            var = self._sources[cmpd]
            var.lb = lb
            var.ub = ub
    else:
        var = self._sources[compound]
        var.lb = lb
        var.ub = ub


def _grb_is_fixed(self, reaction):
    pass

def _grb_export2lp(self, filename):
    filename += ".lp"
    self._model.write(filename)

