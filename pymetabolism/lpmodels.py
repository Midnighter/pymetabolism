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
import copy
import logging

from . import miscellaneous as misc
from .errors import PyMetabolismError


logger = logging.getLogger(__name__)
logger.addHandler(misc.NullHandler())


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
        return super(MetaLPModelFacade, mcls).__new__(mcls, name, bases, attrs)


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
    # deal with Gurobi's annoying log file
    tmp_file = tempfile.mkstemp()[1] # absolute path component
    grb.setParam("LogFile", tmp_file)
    os.remove("gurobi.log")
    # suppress reports to stdout
    grb.setParam("OutputFlag", 0)

    # set class attributes
    attrs["_grb"] = grb
    for key in attrs.iterkeys():
        if key.startswith("_"):
            continue
        try:
            attrs[key] = eval("_grb_" + key)
        except NameError:
            pass

def _grb_initialise(self, name):
    self._model = self._grb.Model(name)
    self._system2grb = dict()
    self._grb2system = dict()
    self._drains = list()
    self._sources = list()


class GurobiFacade(LPModelFacade):
    """
    A unified interface for the Gurobi LP solver.

    Gurobi speaks of variables and linear constraints rather than columns and
    rows. Each variable is represented by a column that is involved with
    multiple linear constraints.
    """

    _gurobipy = None

    def __init__(self, name=""):
        """
        Parameters
        ----------
        name: str (optional)
            Name of the model, largely irrelevant.
        """
        if not self.__class__._gurobipy:
            try:
                self.__class__._gurobipy = __import__(name="gurobipy")
            except ImportError:
                raise ImportError("Gurobi python bindings are required for "\
                        "this functionality, see http://www.gurobi.com/")

        LPModelFacade.__init__(self, name)
        self._gurobipy.setParam("OutputFlag", 0)
        self._model = self._gurobipy.Model(name)
        self._variables = dict()
        self._constraints = dict()

    def __copy__(self):
        cpy = GurobiFacade(self.name)
        cpy._model = self._model.copy()
        cpy._variables = dict()
        for col in cpy._model.getVars():
            cpy._variables[col.getAttr("VarName")] = col
        cpy._constraints = dict()
        for row in cpy._model.getConstrs():
            cpy._constraints[row.getAttr("ConstrName")] = row
        return cpy

    def __deepcopy__(self, memo=dict()):
        return self.__copy__()

    def copy(self):
        return self.__deepcopy__()

    def output(self):
        for (row, constr) in self._constraints.iteritems():
            print row, self._model.getRow(constr)
        print

    def _make_binary(self, names):
        for col in names:
            var = self._variables[col]
            var.setAttr("VType", "B")
        self._model.update()

    def _make_integer(self, names):
        for col in names:
            var = self._variables[col]
            var.setAttr("VType", "I")
        self._model.update()

    def add_column(self, column, coefficients, bounds=tuple()):
        """
        Introduces a new variable.

        Parameters
        ----------
        column: str
            A unique identifier for the column.
        coefficients: iterable
            Iterable over pairs of row identifiers and their coefficients.
        bounds: tuple
            A pair of lower and upper bound on the column variable.

        Notes
        -----
        If the column identifier already exists the bounds are ignored but
        coefficients are not. Look at modify_column_bounds instead.
        """
        var = self._variables.get(column)
        if not var:
            var = self._model.addVar(*bounds, name=column)
            self._variables[column] = var
            self._model.update()
        for (row, factor) in coefficients:
            constraint = self._constraints.get(row)
            if constraint:
                self._model.chgCoeff(constraint, var, factor)
            else:
                self._constraints[row] = self._model.addConstr(
                        self._gurobipy.LinExpr(factor, var),
                        self._gurobipy.GRB.EQUAL, 0.0, name=row)
        self._model.update()

    def add_row(self, row, coefficients):
        """
        Introduces a new constraint.

        Parameters
        ----------
        row: str
            A unique identifier for the row.
        coefficients: iterable
            Iterable over pairs of column identifiers and their coefficients.
        """
        constraint = self._constraints.get(row)
        if not constraint:
            constraint = self._model.addConstr(0.0, self._gurobipy.GRB.EQUAL,
                    0.0, name=row)
            self._constraints[row] = constraint
            self._model.update()
        for (column, factor) in coefficients:
            var = self._variables.get(column)
            if var:
                self._model.chgCoeff(constraint, var, factor)
            else:
                raise PyMetabolismError("modifying coefficient of a"\
                    " non-existant column, please add the column first")
        self._model.update()

    def add_columns(self, columns):
        """
        Bulk introduction of new variables with their boundaries.

        Parameters
        ----------
        columns: iterable
            Iterable over triples of column identifier, lower, and
            upper bound.
        """
        for (column, lb, ub) in columns:
            if not column in self._variables:
                self._variables[column] = self._model.addVar(lb, ub,
                        name=column)
        self._model.update()

    def add_rows(self, rows):
        """
        Bulk introduction of new constraints.

        Parameters
        ----------
        rows: iterable
            Iterable over row identifiers.
        """
        for row in rows:
            if not row in self._constraints:
                self._constraints[row] = self._model.addConstr(0.0,
                        self._gurobipy.GRB.EQUAL, 0.0, name=row.name)
        self._model.update()

    def delete_columns(self, columns):
        """
        Removes column(s) from the model.

        Parameters
        ----------
        columns: iterable or str
            Name or iterable over the column names to be removed.
        """
        if hasattr(columns, "__iter__"):
            for name in columns:
                self._model.remove(self._variables.pop(name))
        else:
            self._model.remove(self._variables.pop(columns))
        self._model.update()

    def delete_row(self, rows):
        """
        Removes row(s) from the model.

        Parameters
        ----------
        rows: iterable or str
            Name or iterable over the row name(s) to be removed.
        """
        if hasattr(rows, "__iter__"):
            for name in rows:
                self._model.remove(self._constraints.pop(name))
        else:
            self._model.remove(self._variables.pop(rows))
        self._model.update()

    def get_columns(self, row=None, coefficients=False):
        """
        Parameters
        ----------
        row: str (optional)
            The name of a single row in the model.
        coefficients: bool (optional)
            In combination with row this will return the respective coefficient
            in each column.

        Returns
        -------
        iterator:
            An iterator over all column names or with the optional parameter
            row, an iterator over all columns with non-zero coefficients the row
            participates in.
        """
        if row:
            lin_expr = self._model.getRow(self._constraints[row])
            if coefficients:
                return ((lin_expr.getVar(i).getAttr("VarName"),
                    lin_expr.getCoeff(i)) for i in xrange(lin_expr.size()))
            else:
                return (lin_expr.getVar(i).getAttr("VarName") for i in
                    xrange(lin_expr.size()))
        else:
            return self._variables.iterkeys()

    def get_rows(self, column=None, coefficients=False):
        """
        Parameters
        ----------
        column: str (optional)
            The name of a single column in the model.
        coefficients: bool (optional)
            In combination with column this will return the respective coefficient
            in each row.

        Returns
        -------
        iterator:
            An iterator over all row names or with the optional parameter
            column, an iterator over all rows with non-zero coefficients the
            column participates in.
        """
        if column:
            grb_column = self._model.getCol(self._variables[column])
            if coefficients:
                return ((grb_column.getConstr(i).getAttr("ConstrName"),
                    grb_column.getCoeff(i)) for i in xrange(grb_column.size()))
            else:
                return (grb_column.getConstr(i).getAttr("ConstrName") for i in
                    xrange(grb_column.size()))
        else:
            return self._constraints.iterkeys()

    def modify_column_coefficients(self, column, coefficients):
        """
        Modify a number of coefficients affecting one column.

        Parameters
        ----------
        column: str
            Name of the column variable to be modified.
        coefficients: iterable
            Iterable over pairs of row identifiers and their coefficients.
        """
        var = self._variables[column]
        for (row, factor) in coefficients:
            self._model.chgCoeff(self._constraints[row], var, factor)
        self._model.update()

    def modify_row_coefficients(self, row, coefficients):
        """
        Modify coefficients affecting a number of columns.

        Parameters
        ----------
        row: str
            Name of the row to be modified.
        coefficients: iterable
            Iterable over pairs of column identifiers and their coefficients.
        """
        constraint = self._constraints[row]
        for (column, factor) in coefficients.iteritems():
            self._model.chgCoeff(constraint, self._variables[column], factor)
        self._model.update()

    def modify_column_bounds(self, columns):
        """
        Modifies the lower and upper bounds of variable(s).

        Parameters
        ----------
        columns: iterable
            Iterable over triples of column identifiers and their lower and
            upper bounds.
        """
        for (name, lb, ub) in columns:
            var = self._variables[name]
            var.setAttr("LB", lb)
            var.setAttr("UB", ub)
        self._model.update()

    def modify_row_bounds(self, rows):
        """
        Modifies the lower and upper bounds of a particular constraint.

        Parameters
        ----------
        rows: iterable
            Iterable over triples of row identifiers and their lower and
            upper bounds.
        """
        raise NotImplementedError("Gurobi does not support bounds for rows")

    def get_column_bounds(self, column):
        """
        Parameters
        ----------
        column: str
            Name of the column whose bounds are to be fetched.

        Returns
        -------
        tuple:
            A pair of the lower and upper bound of the specified column.
        """
        var = self.variables[name]
        return (var.getAttr("LB"), var.getAttr("UB"))

    def get_objective(self, coefficients=False):
        """
        Parameters
        ----------
        coefficients: bool (optional)
            Causes the returned iterator to run over pairs of column name and
            absolute weight in the objective function.

        Returns
        -------
        iterator:
            Current column(s) that are used as objectives in LP.
        """
        lin_expr = self._model.getObjective()
        if coefficients:
            return ((lin_expr.getVar(i).getAttr("VarName"), lin_expr.getCoeff(i))
                    for i in xrange(lin_expr.size()))
        else:
            return (lin_expr.getVar(i).getAttr("VarName") for i in xrange(lin_expr.size()))

    def set_objective(self, columns):
        """
        Determine the variables that are to be maximized or minimized, and their
        relative factors. If one variable needs to be maximized and the other
        minimized in equal weights they can be made objectives with opposing
        sign.

        Parameters
        ----------
        columns: iterable
            Iterable over pairs of column names and coefficients.
        """
        for (name, factor) in columns:
            self._variables[name].setAttr("Obj", factor)
        self._model.update()

    def _status(self):
        """
        Determine the current status of the Gurobi model.
        """
        status = self._model.getAttr("Status")
        if status == self._gurobipy.GRB.LOADED:
            raise PyMetabolismError("optimize before retrieving the objective value")
        elif status == self._gurobipy.GRB.OPTIMAL:
            pass
        elif status == self._gurobipy.GRB.INFEASIBLE:
            raise PyMetabolismError("model is infeasible")
        elif status == self._gurobipy.GRB.INF_OR_UNBD:
            raise PyMetabolismError("model is infeasible or unbounded")
        elif status == self._gurobipy.GRB.UNBOUNDED:
            raise PyMetabolismError("model is unbounded")
        elif status == self._gurobipy.GRB.CUTOFF:
            raise PyMetabolismError("model solution is worse than provided cut-off")
        elif status == self._gurobipy.GRB.ITERATION_LIMIT:
            raise PyMetabolismError("iteration limit exceeded")
        elif status == self._gurobipy.GRB.NODE_LIMIT:
            raise PyMetabolismError("node limit exceeded")
        elif status == self._gurobipy.GRB.TIME_LIMIT:
            raise PyMetabolismError("time limit exceeded")
        elif status == self._gurobipy.GRB.SOLUTION_LIMIT:
            raise PyMetabolismError("solution limit reached")
        elif status == self._gurobipy.GRB.INTERRUPTED:
            raise PyMetabolismError("optimization process was interrupted")
        elif status == self._gurobipy.GRB.SUBOPTIMAL:
            raise PyMetabolismError("solution is suboptimal")
        elif status == self._gurobipy.GRB.NUMERIC:
            raise PyMetabolismError("optimization aborted due to numeric difficulties")

    def get_objective_value(self):
        """
        Returns
        -------
        float:
            Current value of the objective if available.

        Warnings
        --------
        A number of different kinds of exceptions may be raised and inform about
        the status of an available solution.
        """
        self._status()
        try:
            return self._model.getAttr("ObjVal")
        except AttributeError:
            pass

    def get_solution_vector(self, columns=None):
        """
        Parameters
        ----------
        columns: iterable (optional)
            Iterable over column names.

        Returns
        -------
        iterable:
            Iterator over pairs of column identifier and variable value.
        """
        self._status()
        if columns:
            if hasattr(rows, "__iter__"):
                return ((name, self._variables[name].getAttr("X"))\
                        for name in columns)
            else:
                return (columns, self._variables[columns].getAttr("X")) 
        else:
            return ((name, var.getAttr("X")) for (name, var) in\
                    self._variables.iteritems())


    def optimize(self, maximize=True):
        """
        Parameters
        ----------
        maximize: bool (optional)
            Indicates whether the current objective(s) should be maximized or
            minimized.
        """
        if maximize:
            self._model.setAttr("ModelSense", -1)
        else:
            self._model.setAttr("ModelSense", 1)
        self._model.optimize()

    def export2lp(self, filename):
        """
        This was mostly for debugging.
        """
        filename += ".lp"
        self._model.write(filename)

