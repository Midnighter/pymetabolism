#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
====================
LP Solver Interfaces
====================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
    Nikolaus Sonnenschein
:Date:
    2011-03-28
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    lpmodel.py
"""


import warnings
import numpy
try:
    import gurobipy
except ImportError:
    pass
try:
    import glpk
except ImportError:
    pass


class LPModelFacade(object):
    """
    Abstract base class that presents the interface a subclass must implement
    for a specific linear programming solver.
    """

    def __init__(self):
        """
        Warning
        -------
        Instantiation of this class will work (for inheritance) but none of the
        methods are implemented.
        """
        object.__init__(self)

    def add_column(self, name, coefficients):
        """
        Introduces a new variable.

        Parameters
        ----------
        name: str
            Identifier for the column.
        coefficients: dict or iterable
            key-value pairs of row names and their coefficients.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def add_row(self, name, coefficients):
        """
        Introduces a new variable.

        Parameters
        ----------
        name: str
            Identifier for the row.
        coefficients: dict or iterable
            key-value pairs of column names and their coefficients.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def add_columns(self, variables):
        """
        Introduces new variables.

        Parameters
        ----------
        variables: dict or iterable
            dict of dict or an iterable with pairs of variable name and a dict
            of coefficients.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def add_rows(self, variables):
        """
        Introduces new variables.

        Parameters
        ----------
        variables: dict or iterable
            dict of dict or an iterable with pairs of variable name and a dict
            of coefficients.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def delete_column(self, name):
        """
        Removes a variable from the model.

        Parameters
        ----------
        name: str
            Name of the column variable to be removed.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def delete_row(self, name):
        """
        Removes a variable from the model.

        Parameters
        ----------
        name: str
            Name of the row variable to be removed.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def delete_columns(self, variables):
        """
        Removes variables from the model.

        Parameters
        ----------
        variables: iterable
            An iterable that contains the column names as strings.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def delete_rows(self, variables):
        """
        Removes variables from the model.

        Parameters
        ----------
        variables: iterable
            An iterable that contains the row names as strings.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def get_column_names(self):
        """
        Returns
        -------
        list:
            Strings identifying column names.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def get_row_names(self):
        """
        Returns
        -------
        list:
            Strings identifying row names.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def modify_column_bounds(self, name, constraints):
        """
        Modifies the lower and upper bounds of a particular column.

        Parameters
        ----------
        name: str
            Name of the column variable constraints to be modified.
        constraints: dict or iterable:
            Key-value pairs of row names and a pair of lower and upper
            bound.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def modify_row_bounds(self, name, constraints):
        """
        Modifies the lower and upper bounds of a particular row.

        Parameters
        ----------
        name: str
            Name of the row variable constraints to be modified.
        constraints: dict or iterable:
            Key-value pairs of column names and a pair of lower and upper
            bound.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def set_objective(self, variables):
        """
        Modifies the lower and upper bounds of a particular row.

        Parameters
        ----------
        variables: iterable
            Iterable containing column names that should be optimization
            objectives.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def optimize(self, maximize=True):
        """
        Maximizes (or minimizes) the objective(s).
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")


class GurobiLPModelFacade(LPModelFacade):
    """
    A unified interface for the Gurobi LP solver.

    Gurobi speaks of variables and linear constraints rather than rows and
    columns. Each variable is involved in a column with multiple linear
    constraints.
    """

    def __init__(self, name=""):
        """
        Parameters
        ----------
        name: str
            Optional name of the model.
        """
        LPModelFacade.__init__(self)
        self._model = gurobipy.Model(name)
        self._variables = dict()
        self._constraints = dict()

    def _add_variable(self, name, coefficients, bounds):
        """
        """
        if name in self._variables:
            return
        self._variables[name] = self._model.addVar(*bounds, name=name)
        for (row, factor) in coefficients.iteritems():
            if row not in self._constraints:
                self._constraints[row] = self._model.addConstr(gurobipy.LinExpr(factor,
                    self._variables[name]), gurobipy.GRB.EQUAL, 0.0)

    def _add_constraint(self, name, coefficients):
        """
        """
        if name in self._constraints:
            self._model.remove(self._constraints.pop(name))
        self._constraints[name] = self._model.addConstr(gurobipy.LinExpr(
                coefficients.values(), [self._variables[column] for column in
                    coefficients.iterkeys()]), gurobipy.GRB.EQUAL, 0.0)

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
        if hasattr(name, "__iter__"):
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
        if hasattr(name, "__iter__"):
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
        if not isinstance(variables, dict):
            variables = dict(variables)
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
        status = self._model.getAttr("Status")
        if status == gurobipy.GRB.LOADED:
            warnings.warn("optimize before retreiving the objective value")
        elif status == gurobipy.GRB.OPTIMAL:
            return self._model.getAttr("ObjVal")
        elif status == gurobipy.GRB.INFEASIBLE:
            warnings.warn("model is infeasible")
        elif status == gurobipy.GRB.INF_OR_UNBD:
            warnings.warn("model is infeasible or unbounded")
        elif status == gurobipy.GRB.UNBOUNDED:
            warnings.warn("model is unbounded")
        elif status == gurobipy.GRB.CUTOFF:
            warnings.warn("model solution is worse than provided cut-off")
        elif status == gurobipy.GRB.ITERATION_LIMIT:
            warnings.warn("iteration limit exceeded")
        elif status == gurobipy.GRB.NODE_LIMIT:
            warnings.warn("node limit exceeded")
        elif status == gurobipy.GRB.TIME_LIMIT:
            warnings.warn("time limit exceeded")
        elif status == gurobipy.GRB.SOLUTION_LIMIT:
            return self._model.getAttr("ObjVal")
        elif status == gurobipy.GRB.INTERRUPTED:
            warnings.warn("optimization process was interrupted")
        elif status == gurobipy.GRB.SUBOPTIMAL:
            warnings.warn("solution is suboptimal")
            return self._model.getAttr("ObjVal")
        elif status == gurobipy.GRB.NUMERIC:
            warnings.warn("optimization aborted due to numeric difficulties")

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


