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


import copy
import logging

from . import miscellaneous as misc


logger = logging.getLogger(__name__)
logger.addHandler(misc.NullHandler())


class LPModelFacade(object):
    """
    Abstract base class that presents the interface a subclass must implement
    for a specific linear programming solver.
    """

    def __init__(self, name=""):
        """
        Warning
        -------
        Instantiation of this class will work (for inheritance) but none of the
        methods are implemented.
        """
        object.__init__(self)
        self.name = name

    def __str__(self):
        return self.name

    def add_column(self, name, coefficients, bounds=tuple()):
        """
        Introduces a new variable.

        Parameters
        ----------
        name: str
            Identifier for the column.
        coefficients: dict
            key-value pairs of row names and their coefficients.
        """

    def add_row(self, name, coefficients):
        """
        Introduces a new constraint.

        Parameters
        ----------
        name: str
            Identifier for the row.
        coefficients: dict
            key-value pairs of column names and their coefficients.
        """

    def add_columns(self, variables):
        """
        Introduces new variables.

        Parameters
        ----------
        variables: iterable
            iterable with triples of variable name, a dict of coefficients, and
            a tuple with lower and upper bound.
        """

    def add_rows(self, variables):
        """
        Introduces new constraints.

        Parameters
        ----------
        variables: dict
            dict of dict that maps constraint name(s) to a dict of coefficients.
        """

    def delete_column(self, column):
        """
        Removes a column from the model.

        Parameters
        ----------
        column: iterable or str
            Name or iterable over the column name(s) to be removed.
        """

    def delete_row(self, row):
        """
        Removes a row from the model.

        Parameters
        ----------
        row: iterable or str
            Name or iterable over the row name(s) to be removed.
        """

    def get_column_names(self, row=None, coefficients=False):
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

    def get_row_names(self, column=None, coefficients=False):
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

    def modify_column_coefficients(self, name, coefficients):
        """
        Modify a number of coefficients affecting one column.

        Parameters
        ----------
        name: str
            Name of the column variable to be modified.
        constraints: dict
            Map from row names to a pair of lower and upper
            bounds.
        """
        var = self._variables[name]
        for (row, factor) in constraints.iteritems():
            self._model.chgCoeff(self._constraints[row], var, factor)
        self._model.update()

    def modify_row_coefficients(self, name, coefficients):
        """
        Modify coefficients affecting a number of columns.

        Parameters
        ----------
        name: str
            Name of the row to be modified.
        constraints: dict or iterable
            Map of column name(s) to the new coefficient.
        """
        constraint = self._constraints[name]
        for (column, factor) in constraints.iteritems():
            self._model.chgCoeff(constraint, self._variables[column], factor)
        self._model.update()

    def modify_column_bounds(self, variables):
        """
        Modifies the lower and upper bounds of variable(s).

        Parameters
        ----------
        variables: dict
            Map of column name to a pair with lower and upper bound.
        """
        for (name, bounds) in variables.iteritems():
            self._variables[name].setAttr("LB", bounds[0])
            self._variables[name].setAttr("UB", bounds[1])
        self._model.update()

    def modify_row_bounds(self, name, bounds):
        """
        Modifies the lower and upper bounds of a particular row.

        Parameters
        ----------
        name: str
            Name of the row variable constraints to be modified.
        constraints: dict or iterable
            Key-value pairs of column names and a pair of lower and upper
            bound.
        """
        raise NotImplementedError("Gurobi does not support bounds for rows")
    
    def get_column_bounds(self, name):
        return (self._variables[name].getAttr("LB"),
                self._variables[name].getAttr("UB"))

    def get_objective(self, coefficient=False):
        """
        Parameters
        ----------
        coefficient: bool (optional)
            Causes the returned iterator to run over pairs of column name and
            absolute weight in the objective function.

        Returns
        -------
        iterator:
            Current column(s) that are used as objectives in LP.
        """
        lin_expr = self._model.getObjective()
        if coefficient:
            return ((lin_expr(i).getAttr("VarName"), lin_expr(i).getCoeff())
                    for i in xrange(lin_expr.size()))
        else:
            return (lin_expr(i).getAttr("VarName") for i in xrange(lin_expr.size()))

    def set_objective(self, variables):
        """
        Determine the variables that are to be maximized or minimized, and their
        relative factors. If one variable needs to be maximized and the other
        minimized in equal weights they can be made objectives with opposing
        sign.

        Parameters
        ----------
        variables: dict
            Map from column name(s) to factor(s).
        """
        for (name, factor) in variables.iteritems():
            self._variables[name].setAttr("Obj", factor)
        self._model.update()

    def get_objective_value(self):
        """
        Returns
        -------
        float:
            Current value of the objective if available.

        Warnings
        --------
        A number of different kinds of warnings may be issued and inform about
        the status of an available solution.
        """

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


        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")


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

    def _add_variable(self, name, coefficients, bounds):
        """
        """
        if name in self._variables:
            return
        self._variables[name] = self._model.addVar(*bounds, name=name)
        self._model.update()
        for (row, factor) in coefficients.iteritems():
            if row not in self._constraints:
                self._constraints[row] = self._model.addConstr(
                        self._gurobipy.LinExpr(factor, self._variables[name]),
                        self._gurobipy.GRB.EQUAL, 0.0, name=row)
                self._model.update()

    def _add_constraint(self, name, coefficients):
        """
        """
        if name in self._constraints:
            lin_expr = self._model.getRow(self._constraints[name])
            msg = [(lin_expr.getCoeff(i), lin_expr.getVar(i).getAttr("VarName"))
                    for i in xrange(lin_expr.size())]
            msg = " ".join(str(c) for expr in msg for c in expr)
            msg = "replacing existing row '%s'\n%s" % (name, msg)
            logger.warn(msg)
            self._model.remove(self._constraints.pop(name))
        self._constraints[name] = self._model.addConstr(self._gurobipy.LinExpr(
                coefficients.values(), [self._variables[column] for column in
                coefficients.iterkeys()]), self._gurobipy.GRB.EQUAL, 0.0,
                name=name)
        self._model.update()

    def add_column(self, name, coefficients, bounds=tuple()):
        """
        Introduces a new variable.

        Parameters
        ----------
        name: str
            Identifier for the column.
        coefficients: dict
            key-value pairs of row names and their coefficients.
        """
        self._add_variable(name, coefficients, bounds)
        var = self._variables[name]
        for (row, factor) in coefficients.iteritems():
            self._model.chgCoeff(self._constraints[row], var, factor)
        self._model.update()

    def add_row(self, name, coefficients):
        """
        Introduces a new constraint.

        Parameters
        ----------
        name: str
            Identifier for the row.
        coefficients: dict
            key-value pairs of column names and their coefficients.
        """
        self._add_constraint(name, coefficients)
        self._model.update()

    def add_columns(self, variables):
        """
        Introduces new variables.

        Parameters
        ----------
        variables: iterable
            iterable with triples of variable name, a dict of coefficients, and
            a tuple with lower and upper bound.
        """
        add_func = self._add_variable
        chg_func = self._model.chgCoeff
        for (name, coefficients, bounds) in variables:
            add_func(name, coefficients, bounds)
            var = self._variables[name]
            for (row, factor) in coefficients.iteritems():
                chg_func(self._constraints[row], var, factor)
        self._model.update()

    def add_rows(self, variables):
        """
        Introduces new constraints.

        Parameters
        ----------
        variables: dict
            dict of dict that maps constraint name(s) to a dict of coefficients.
        """
        add_func = self._add_constraint
        for (name, coefficients) in variables.iteritems():
            add_func(name, coefficients)
        self._model.update()

    def delete_column(self, column):
        """
        Removes a column from the model.

        Parameters
        ----------
        column: iterable or str
            Name or iterable over the column name(s) to be removed.
        """
        if hasattr(column, "__iter__"):
            for name in column:
                self._model.remove(self._variables.pop(name))
        else:
            self._model.remove(self._variables.pop(column))
        self._model.update()

    def delete_row(self, row):
        """
        Removes a row from the model.

        Parameters
        ----------
        row: iterable or str
            Name or iterable over the row name(s) to be removed.
        """
        if hasattr(row, "__iter__"):
            for name in row:
                self._model.remove(self._constraints.pop(name))
        else:
            self._model.remove(self._variables.pop(row))
        self._model.update()

    def get_column_names(self, row=None, coefficients=False):
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

    def get_row_names(self, column=None, coefficients=False):
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

    def modify_column_coefficients(self, name, coefficients):
        """
        Modify a number of coefficients affecting one column.

        Parameters
        ----------
        name: str
            Name of the column variable to be modified.
        constraints: dict
            Map from row names to a pair of lower and upper
            bounds.
        """
        var = self._variables[name]
        for (row, factor) in constraints.iteritems():
            self._model.chgCoeff(self._constraints[row], var, factor)
        self._model.update()

    def modify_row_coefficients(self, name, coefficients):
        """
        Modify coefficients affecting a number of columns.

        Parameters
        ----------
        name: str
            Name of the row to be modified.
        constraints: dict or iterable
            Map of column name(s) to the new coefficient.
        """
        constraint = self._constraints[name]
        for (column, factor) in constraints.iteritems():
            self._model.chgCoeff(constraint, self._variables[column], factor)
        self._model.update()

    def modify_column_bounds(self, variables):
        """
        Modifies the lower and upper bounds of variable(s).

        Parameters
        ----------
        variables: dict
            Map of column name to a pair with lower and upper bound.
        """
        for (name, bounds) in variables.iteritems():
            self._variables[name].setAttr("LB", bounds[0])
            self._variables[name].setAttr("UB", bounds[1])
        self._model.update()

    def modify_row_bounds(self, name, bounds):
        """
        Modifies the lower and upper bounds of a particular row.

        Parameters
        ----------
        name: str
            Name of the row variable constraints to be modified.
        constraints: dict or iterable
            Key-value pairs of column names and a pair of lower and upper
            bound.
        """
        raise NotImplementedError("Gurobi does not support bounds for rows")

    def get_column_bounds(self, name):
        """
        Parameters
        ----------
        name: str
            Name of the column whose bounds are to be fetched.

        Returns
        -------
        tuple:
            A pair of the lower and upper bound of the specified column.
        """
        return (self._variables[name].getAttr("LB"),
                self._variables[name].getAttr("UB"))

    def get_objective(self, coefficient=False):
        """
        Parameters
        ----------
        coefficient: bool (optional)
            Causes the returned iterator to run over pairs of column name and
            absolute weight in the objective function.

        Returns
        -------
        iterator:
            Current column(s) that are used as objectives in LP.
        """
        lin_expr = self._model.getObjective()
        if coefficient:
            return ((lin_expr.getVar(i).getAttr("VarName"), lin_expr.getCoeff(i))
                    for i in xrange(lin_expr.size()))
        else:
            return (lin_expr.getVar(i).getAttr("VarName") for i in xrange(lin_expr.size()))

    def set_objective(self, variables):
        """
        Determine the variables that are to be maximized or minimized, and their
        relative factors. If one variable needs to be maximized and the other
        minimized in equal weights they can be made objectives with opposing
        sign.

        Parameters
        ----------
        variables: dict
            Map from column name(s) to factor(s).
        """
        for (name, factor) in variables.iteritems():
            self._variables[name].setAttr("Obj", factor)
        self._model.update()

    def _status(self):
        status = self._model.getAttr("Status")
        if status == self._gurobipy.GRB.LOADED:
            logger.warn("optimize before retrieving the objective value")
        elif status == self._gurobipy.GRB.OPTIMAL:
            return self._model.getAttr("ObjVal")
        elif status == self._gurobipy.GRB.INFEASIBLE:
            logger.warn("model is infeasible")
        elif status == self._gurobipy.GRB.INF_OR_UNBD:
            logger.warn("model is infeasible or unbounded")
        elif status == self._gurobipy.GRB.UNBOUNDED:
            logger.warn("model is unbounded")
        elif status == self._gurobipy.GRB.CUTOFF:
            logger.warn("model solution is worse than provided cut-off")
        elif status == self._gurobipy.GRB.ITERATION_LIMIT:
            logger.warn("iteration limit exceeded")
        elif status == self._gurobipy.GRB.NODE_LIMIT:
            logger.warn("node limit exceeded")
        elif status == self._gurobipy.GRB.TIME_LIMIT:
            logger.warn("time limit exceeded")
        elif status == self._gurobipy.GRB.SOLUTION_LIMIT:
            return self._model.getAttr("ObjVal")
        elif status == self._gurobipy.GRB.INTERRUPTED:
            logger.warn("optimization process was interrupted")
        elif status == self._gurobipy.GRB.SUBOPTIMAL:
            logger.warn("solution is suboptimal")
            return self._model.getAttr("ObjVal")
        elif status == self._gurobipy.GRB.NUMERIC:
            logger.warn("optimization aborted due to numeric difficulties")
        return status

    def get_objective_value(self):
        """
        Returns
        -------
        float:
            Current value of the objective if available.

        Warnings
        --------
        A number of different kinds of warnings may be issued and inform about
        the status of an available solution.
        """
        self._status()
        try:
            return self._model.getAttr("ObjVal")
        except AttributeError:
            pass

    def get_solution_vector(self):
        self._status()
        return ((column, var.getAttr("X")) for (column, var) in
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
        return self._status()

    def export2lp(self, filename):
        filename += ".lp"
        self._model.write(filename)

