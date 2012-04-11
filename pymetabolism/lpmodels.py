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
    grb.setParam("Threads", options.num_proc)
    # set the feasability tolerance (smaller is more accurate but harder)
    grb.setParam("FeasibilityTol", options.numeric_threshold)

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
    attrs["__init__"] = _grb___init__
#    attrs["__copy__"] = _grb___copy__
#    attrs["__deepcopy__"] = _grb___deepcopy__
    attrs["_add_compound"] = _grb__add_compound
    attrs["_change_participation"] = _grb__change_participation
    attrs["_make_binary"] = _grb__make_binary
    attrs["_make_integer"] = _grb__make_integer
    attrs["_add_reaction"] = _grb__add_reaction
    attrs["_change_coefficients"] = _grb__change_coefficients
    attrs["_adjust_bounds"] = _grb__adjust_bounds
    attrs["_bounds"] = _grb__bounds
    attrs["_fixed"] = _grb__fixed
    attrs["_add_transport"] = _grb__add_transport
    attrs["_add_source"] = _grb__add_source
    attrs["_add_drain"] = _grb__add_drain
    attrs["_var2reaction"] = _grb__var2reaction
    attrs["_reset_objective"] = _grb__reset_objective
    attrs["_status"] = _grb__status
    attrs["_flux"] = _grb__flux
    attrs["_reduced_cost"] = _grb__reduced_cost

def _grb___init__(self, name):
    self._model = self._grb.Model(name)
    self._rxn2var = dict()
    self._var2rxn = dict()
    self._rev2var = dict()
    self._var2rev = dict()
    self._cmpd2cnstrnt = dict()
    self._cnstrnt2cmpd = dict()
    self._sources = dict()
    self._drains = dict()
    self._objective = dict()
    self._tmp_lb = dict()

#def _grb___copy__(self):
#    # TODO
#    cpy = self.__class__(self.name)
#    cpy._model = self._model.copy()
#    cpy._system2grb = dict()
#    for col in cpy._model.getVars():
#        cpy._system2grb[col.getAttr("VarName")] = col
#    cpy._system2grb = dict()
#    for row in cpy._model.getConstrs():
#        cpy._system2grb[row.getAttr("ConstrName")] = row
#    return cpy
#
#def _grb___deepcopy__(self, memo=dict()):
#    # TODO
#    return self.__copy__()
#
#def _grb_copy(self):
#    # TODO
#    return self.__deepcopy__()

def _grb__add_compound(self, compound):
    if compound in self._cmpd2cnstrnt:
        return False
    cnstrnt = self._model.addConstr(0.0, self._grb.GRB.EQUAL,
            0.0, name=str(compound))
    self._cmpd2cnstrnt[compound] = cnstrnt
    self._cnstrnt2cmpd[cnstrnt] = compound
    return True

def _grb__change_participation(self, compound, coefficients):
    cnstrnt = self._cmpd2cnstrnt[compound]
    for (rxn, factor) in coefficients:
        # will throw KeyError if reaction doesn't exist yet
        var = self._rxn2var[rxn]
        self._model.chgCoeff(cnstrnt, var, factor)
        if rxn.reversible:
            var = self._rev2var[rxn]
            self._model.chgCoeff(cnstrnt, var, -factor)

def _grb_add_compound(self, compound, coefficients=None):
    if hasattr(compound, "__iter__"):
    # we really add multiple compounds
        if hasattr(coefficients, "__iter__"):
            coefficients = itertools.repeat(coefficients)
        changes = [self._add_compound(cmpd) for cmpd in compound]
        if any(changes):
            self._model.update()
        if coefficients is None:
            return
        for (cmpd, coeff_iter) in itertools.izip(compound, coefficients):
            if self._model.getRow(self._cmpd2cnstrnt[cmpd]).size() > 0:
                # compound participates in existing reactions thus was added before
                continue
            self._change_participation(cmpd, coeff_iter)
    else:
        if self._add_compound(compound):
            self._model.update()
            if coefficients is None:
                return
            self._change_participation(compound, coefficients)

def _grb_iter_compounds(self, reaction=None, coefficients=False):
    # updating is the only way currently to return newly added information
    self._model.update()
    if reaction is None:
        return self._cmpd2cnstrnt.iterkeys()
    column = self._model.getCol(self._rxn2var[reaction])
    if coefficients:
        return ((self._cnstrnt2cmpd[column.getConstr(i)], column.getCoeff(i))\
                for i in range(column.size()))
    else:
        return (self._cnstrnt2cmpd[column.getConstr(i)]\
                for i in range(column.size()))

def _grb_modify_compound_coefficients(self, compound, coefficients):
    if hasattr(compound, "__iter__"):
        if hasattr(coefficients, "__iter__"):
            coefficients = itertools.repeat(coefficients)
        for (cmpd, coeff_iter) in itertools.izip(compound, coefficients):
            self._change_participation(cmpd, coeff_iter)
    else:
        self._change_participation(compound, coefficients)

def _grb_free_compound(self, compound):
    pass

def _grb_knockout_compound(self, compound):
    cnstrnt = self._cmpd2cnstrnt[compound]
    lin_expr = self._model.getRow(cnstrnt)
    for i in range(lin_expr.size()):
        var = lin_expr.getVar(i)
        var.lb = 0.0
        var.ub = 0.0

def _grb__del_compound(self, compound):
    cnstrnt = self._cmpd2cnstrnt.pop(compound)
    del self._cnstrnt2cmpd[cnstrnt]
    self._model.remove(cnstrnt)

def _grb_delete_compound(self, compound):
    if hasattr(compound, "__iter__"):
        for cmpd in compound:
            self._del_compound(cmpd)
    else:
        self._del_compound(compound)
    self._model.update()

def _grb__make_binary(self, reaction):
    if hasattr(reaction, "__iter__"):
        for rxn in reaction:
            var = self._rxn2var[rxn]
            var.vType = "B"
            if rxn.reversible:
                var = self._rev2var[rxn]
                var.vType = "B"
    else:
        var = self._rxn2var[reaction]
        var.vType = "B"
        if reaction.reversible:
            var = self._rev2var[reaction]
            var.vType = "B"

def _grb__make_integer(self, reaction):
    if hasattr(reaction, "__iter__"):
        for rxn in reaction:
            var = self._rxn2var[rxn]
            var.vType = "I"
            if rxn.reversible:
                var = self._rev2var[rxn]
                var.vType = "I"
    else:
        var = self._rxn2var[reaction]
        var.vType = "I"
        if reaction.reversible:
            var = self._rev2var[reaction]
            var.vType = "I"

def _grb__add_reaction(self, reaction, lb, ub):
    if self._rxn2var.has_key(reaction):
        return False
    if reaction.reversible:
        # we rely on lb being numeric here due to default options
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
    return True

def _grb__change_coefficients(self, reaction, coefficients):
    var = self._rxn2var[reaction]
    for (cmpd, factor) in coefficients:
        cnstrnt = self._cmpd2cnstrnt[cmpd]
        self._model.chgCoeff(cnstrnt, var, factor)
    if reaction.reversible:
        var = self._rev2var[reaction]
        for (cmpd, factor) in coefficients:
            cnstrnt = self._cmpd2cnstrnt[cmpd]
            self._model.chgCoeff(cnstrnt, var, -factor)

def _grb_add_reaction(self, reaction, coefficients=None, lb=None, ub=None):
    if lb is None:
        lb = options.lower_bound
    if ub is None:
        ub = options.upper_bound
    if hasattr(reaction, "__iter__"):
    # we really add multiple reactions
        if hasattr(lb, "__iter__"):
            lb_iter = lb
        else:
            lb_iter = itertools.repeat(lb)
        if hasattr(ub, "__iter__"):
            ub_iter = ub
        else:
            ub_iter = itertools.repeat(ub)
        changes = [self._add_reaction(rxn, lb, ub) for (rxn, lb, ub)\
                in itertools.izip(reaction, lb_iter, ub_iter)]
        if any(changes):
            self._model.update()
        if coefficients is None:
            return
        # need to find out if we are dealing with a nested list or not
        if not (isinstance(coefficients, list) and isinstance(coefficients[0],
                list)):
            coefficients = itertools.repeat(coefficients)
        for (rxn, coeff_iter) in itertools.izip(reaction, coefficients):
            changes = [self._add_compound(pair[0]) for pair in coeff_iter]
            if any(changes):
                self._model.update()
            if self._model.getCol(self._rxn2var[rxn]).size() > 0:
                # reaction has constraints and was added before
                continue
            self._change_coefficients(rxn, coeff_iter)
    else:
        if self._add_reaction(reaction, lb, ub):
            self._model.update()
            if coefficients is None:
                return
            changes = [self._add_compound(pair[0]) for pair in coeff_iter]
            if any(changes):
                self._model.update()
            self._change_coefficients(rxn, coeff_iter)

def _grb_iter_reactions(self, compound=None, coefficients=False):
    # updating is the only way currently to return newly added information
    self._model.update()
    if not compound:
        return self._rxn2var.iterkeys()
    lin_expr = self._model.getRow(self._cmpd2cnstrnt[compound])
    if coefficients:
        return ((self._var2reaction(lin_expr.getVar(i)), lin_expr.getCoeff(i))\
                for i in range(lin_expr.size()))
    else:
        return (self._var2reaction(lin_expr.getVar(i))\
                for i in range(lin_expr.size()))

def _grb_modify_reaction_coefficients(self, reaction, coefficients):
    # we allow for lazy updating of the model here (better not be a bug)
    if hasattr(reaction, "__iter__"):
        if hasattr(coefficients, "__iter__"):
            coefficients = itertools.repeat(coefficients)
        for (rxn, coeff_iter) in itertools.izip(reaction, coefficients):
            self._change_coefficients(rxn, coeff_iter)
    else:
        self._change_coefficients(reaction, coefficients)

def _grb__adjust_bounds(self, reaction, lb, ub):
    """
    Adjust the lower and upper bound for a reaction.

    Reversible reactions are treated specially since bounds may be split for
    both directions.
    """
    numeric_ub = not ub is None
    numeric_lb = not lb is None
    if numeric_ub and numeric_lb and ub < lb:
        raise PyMetabolismError("Trying to set an upper bound that is smaller"\
        " than the lower bound for '%s'.", str(reaction))
    var = self._rxn2var[reaction]
    if reaction.reversible:
        var_rev = self._rev2var[reaction]
        if numeric_ub:
            var.ub = ub
            var_rev.ub = ub
        if numeric_lb:
            if lb < 0.0:
                var_rev.lb = 0.0
                var_rev.ub = abs(lb)
                var.lb = 0.0
            else:
                var_rev.lb = lb
                var.lb = lb
    else:
        if numeric_ub:
            var.ub = ub
        if numeric_lb:
            var.lb = lb

def _grb_modify_reaction_bounds(self, reaction, lb=None, ub=None):
    # we allow for lazy updating of the model here (better not be a bug)
    if hasattr(reaction, "__iter__"):
        # we really modify multiple reactions
        if hasattr(lb, "__iter__"):
            lb_iter = lb
        else:
            lb_iter = itertools.repeat(lb)
        if hasattr(ub, "__iter__"):
            ub_iter = ub
        else:
            ub_iter = itertools.repeat(ub)
        for (rxn, lb, ub) in itertools.izip(reaction, lb_iter, ub_iter):
            self._adjust_bounds(rxn, lb, ub)
    else:
        self._adjust_bounds(reaction, lb, ub)

def _grb__bounds(self, reaction):
    var = self._rxn2var[reaction]
    if reaction.reversible:
        var_rev = self._rev2var[reaction]
        return (-var_rev.ub, var.ub)
    else:
        return (var.lb, var.ub)

def _grb_iter_reaction_bounds(self, reaction=None):
    # updating is the only way currently to return newly added information
    self._model.update()
    # we rely on reversible reactions being treated in unison
    if reaction is None:
        reaction = self._rxn2var.iterkeys()
        return ((rxn, self._bounds(rxn)) for rxn in reaction)
    elif hasattr(reaction, "__iter__"):
        # we really get multiple reactions
        return (self._bounds(rxn) for rxn in reaction)
    else:
        return self._bounds(reaction)

def _grb__fixed(self, reaction):
    var = self.rxn2var[reaction]
    fixed = var.lb == var.ub
    if reaction.reversible:
        var = self.rev2var[reaction]
        fixed &= var.lb == var.ub
    return fixed

def _grb_is_fixed(self, reaction=None):
    # updating is the only way currently to return newly added information
    self._model.update()
    if reaction is None:
        reaction = self._rxn2var.iterkeys()
        return all(self._fixed(rxn) for rxn in reaction)
    elif hasattr(reaction, "__iter__"):
        # we really get multiple reactions
        return all(self._fixed(rxn) for rxn in reaction)
    else:
        return self._bounds(reaction)

def _grb_free_reaction(self, reaction):
    self.modify_reaction_bounds(reaction, lb=-self._grb.GRB.INFINITY,
            ub=self._grb.GRB.INFINITY)

def _grb__del_reaction(self, reaction):
    var = self._rxn2var.pop(reaction)
    del self._var2rxn[var]
    self._model.remove(var)
    if reaction.reversible:
        var = self._rev2var.pop(reaction)
        del self._var2rev[var]
        self._model.remove(var)

def _grb_delete_reaction(self, reaction):
    if hasattr(reaction, "__iter__"):
        for rxn in reaction:
            self._del_reaction(rxn)
    else:
        self._del_reaction(reaction)
    self._model.update()

def _grb__add_transport(self, compound, var, factor):
    if self._model.getCol(var).size() > 0:
        # transport already added
        return
    cnstrnt = self._cmpd2cnstrnt[compound]
    self._model.chgCoeff(cnstrnt, var, factor)

def _grb__add_source(self, compound, lb, ub):
    if compound in self._sources:
        return False
    self._sources[compound] = self._model.addVar(lb, ub,
            name=str(compound) + "_Source")
    return True

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
        changes = [self._add_source(cmpd, lb, ub) for (cmpd, lb, ub)\
                in itertools.izip(compound, lb_iter, ub_iter)]
        if any(changes):
            self._model.update()
        # we allow for lazy updating of the model here (better not be a bug)
        for cmpd in compound:
            var = self._sources[cmpd]
            self._add_transport(cmpd, var, 1.0)
    else:
        if self._add_source(compound, lb, ub):
            self._model.update()
            var = self._sources[compound]
            self._add_transport(compound, var, 1.0)

def _grb__add_drain(self, compound, lb, ub):
    if compound in self._drains:
        return False
    self._drains[compound] = self._model.addVar(lb, ub,
            name=str(compound) + "_Drain")
    return True

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
        changes = [self._add_drain(cmpd, lb, ub) for (cmpd, lb, ub)\
                in itertools.izip(compound, lb_iter, ub_iter)]
        if any(changes):
            self._model.update()
        # we allow for lazy updating of the model here (better not be a bug)
        for cmpd in compound:
            var = self._drains[cmpd]
            self._add_transport(cmpd, var, -1.0)
    else:
        if self._add_drain(compound, lb, ub):
            self._model.update()
            var = self._drains[compound]
            self._add_transport(compound, var, -1.0)

def _grb_set_objective_reaction(self, reaction, factor):
    # we allow for lazy updating of the model here (better not be a bug)
    self._objective = dict()
    if hasattr(reaction, "__iter__"):
        if hasattr(factor, "__iter__"):
            fctr_iter = factor
        else:
            fctr_iter = itertools.repeat(factor)
        for (rxn, factor) in itertools.izip(reaction, fctr_iter):
            self._objective[rxn] = factor
    else:
        self._objective[reaction] = factor

def _grb__var2reaction(self, var):
    return self._var2rxn[var] if var in self._var2rxn else self._var2rev[var]

def _grb_iter_objective_reaction(self, coefficients=False):
    if coefficients:
        return self._objective.iteritems()
    else:
        return self._objective.iterkeys()

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

def _grb__reset_objective(self):
    lin_expr = self._model.getObjective()
    for i in range(lin_expr.size()):
        var = lin_expr.getVar(i).obj = 0.0
    for (rxn, factor) in self._objective.iteritems():
        var = self._rxn2var[rxn]
        var.obj = factor
        var.lb = self._tmp_lb.get(var, var.lb)
        if rxn.reversible:
            var = self._rev2var[rxn]
            var.obj = factor
            var.lb = self._tmp_lb.get(var, var.lb)

def _grb_fba(self, maximize=True):
    self._reset_objective()
    if maximize:
        self._model.modelSense = self._grb.GRB.MAXIMIZE
    else:
        self._model.modelSense = self._grb.GRB.MINIMIZE
    self._model.optimize()

def _grb_parsimonious_fba(self):
    self._reset_objective()
    self._model.modelSense = self._grb.GRB.MAXIMIZE
    self._model.optimize()
    # _status should catch all problems (monitor this)
    self._status()
    lin_expr = self._model.getObjective()
    objective = set([lin_expr.getVar(i) for i in range(lin_expr.size())])
    for (i, var) in enumerate(objective):
        self._tmp_lb[var] = var.lb
        var.lb = var.x
        var.obj = 0.0
    for var in itertools.chain(self._rxn2var.itervalues(),
            self._rev2var.itervalues()):
        if not var in objective:
            var.obj = 1.0
    self._model.modelSense = self._grb.GRB.MINIMIZE
    self._model.optimize()

def _grb__status(self):
    """
    Determine the current status of the Gurobi model.
    """
    status = self._model.status
    if status == self._grb.GRB.LOADED:
        raise PyMetabolismError("optimize before accessing flux information", errorno=status)
    elif status == self._grb.GRB.OPTIMAL:
        pass
    elif status == self._grb.GRB.INFEASIBLE:
        raise PyMetabolismError("model is infeasible", errorno=status)
    elif status == self._grb.GRB.INF_OR_UNBD:
        raise PyMetabolismError("model is infeasible or unbounded", errorno=status)
    elif status == self._grb.GRB.UNBOUNDED:
        raise PyMetabolismError("model is unbounded", errorno=status)
    elif status == self._grb.GRB.CUTOFF:
        raise PyMetabolismError("model solution is worse than provided cut-off", errorno=status)
    elif status == self._grb.GRB.ITERATION_LIMIT:
        raise PyMetabolismError("iteration limit exceeded", errorno=status)
    elif status == self._grb.GRB.NODE_LIMIT:
        raise PyMetabolismError("node limit exceeded", errorno=status)
    elif status == self._grb.GRB.TIME_LIMIT:
        raise PyMetabolismError("time limit exceeded", errorno=status)
    elif status == self._grb.GRB.SOLUTION_LIMIT:
        raise PyMetabolismError("solution limit reached", errorno=status)
    elif status == self._grb.GRB.INTERRUPTED:
        raise PyMetabolismError("optimization process was interrupted", errorno=status)
    elif status == self._grb.GRB.SUBOPTIMAL:
        raise PyMetabolismError("solution is suboptimal", errorno=status)
    elif status == self._grb.GRB.NUMERIC:
        raise PyMetabolismError("optimization aborted due to numeric difficulties", errorno=status)

def _grb_get_objective_value(self, threshold=None):
    # _status should catch all problems (monitor this)
    self._status()
    if threshold is None:
        threshold = options.numeric_threshold
    return sum(self._flux(rxn, threshold) * factor for (rxn, factor)\
            in self._objective.iteritems())

def _grb__flux(self, reaction, threshold):
    flux = self._rxn2var[reaction].x
    if reaction.reversible:
        flux -= self._rev2var[reaction].x
    return flux if abs(flux) > threshold else 0.0

def _grb_iter_flux(self, reaction=None, threshold=None):
    # _status should catch all problems (monitor this)
    self._status()
    if threshold is None:
        threshold = options.numeric_threshold
    if reaction is None:
        return ((rxn, self._flux(rxn, threshold)) for rxn in\
                self._rxn2var.iterkeys())
    elif hasattr(reaction, "__iter__"):
        return (self._flux(rxn, threshold) for rxn in reaction)
    else:
        return self._flux(reaction, threshold)

def _grb__reduced_cost(self, reaction, threshold):
    cost = self._rxn2var[reaction].rc
    if reaction.reversible:
        cost -= self._rev2var[reaction].rc
    return cost if abs(cost) > threshold else 0.0

def _grb_iter_reduced_cost(self, reaction=None, threshold=None):
    # _status should catch all problems (monitor this)
    self._status()
    if threshold is None:
        threshold = options.numeric_threshold
    if reaction is None:
        return ((rxn, self._reduced_cost(rxn, threshold)) for rxn in\
                self._rxn2var.iterkeys())
    elif hasattr(reaction, "__iter__"):
        return (self._reduced_cost(rxn, threshold) for rxn in reaction)
    else:
        return self._reduced_cost(reaction, threshold)

def _grb_iter_shadow_price(self, compound=None):
    # _status should catch all problems (monitor this)
    self._status()
    if compound is None:
        compound = self._rxn2var.iterkeys()
        return ((cmpd, cnstrnt.pi) for (cmpd, cnstrnt) in\
                self._cmpd2cnstrnt.iteritems())
    elif hasattr(compound, "__iter__"):
        return (self._cmpd2cnstrnt[cmpd].pi for cmpd in compound)
    else:
        return self._cmpd2cnstrnt[cmpd].pi

def _grb_export2lp(self, filename):
    filename += ".lp"
    self._model.write(filename)

