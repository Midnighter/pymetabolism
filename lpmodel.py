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


import numpy
try:
    import gurobipy
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

    def _get_variable(self, name):
        """
        """
        return self._variables[name]

    def _get_constraint(self, name):
        """
        """
        return self._constraints[name]

    def _add_variable(self, name, coefficients, bounds):
        """
        """
        if name in self._variables:
            return
        self._variables[name] = self._model.addVar(*bounds, name=name)
        for (row, factor) in coefficients.iteritems():
            if row not in self._rows:
                self._rows[row] = self._model.addConstr(gurobipy.LinExpr(factor,
                    self._columns[name]), gurobipy.GRB.EQUAL, 0.0)

    def _add_constraint(self, name, coefficients):
        """
        """
        if name in self._constraints:
            self._model.remove(self._constraints.pop(name))
        self._constraints[name] = self._model.addConstr(gurobipy.LinExpr(
                coefficients.values(), map(self._get_variable,
                coefficients.keys())), gurobipy.GRB.EQUAL, 0.0)

    def add_column(self, name, coefficients, bounds=tuple()):
        """
        Introduces a new variable.

        Parameters
        ----------
        name: str
            Identifier for the column.
        coefficients: dict or iterable
            key-value pairs of row names and their coefficients.
        """
        if not isinstance(coefficients, dict):
            coefficients = dict(coefficients)
        self._add_variable(name, coefficients, bounds)
        grb_column = self._columns[name]
        grb_column.addTerms(coefficients.values(), map(self._get_constraint,
                coefficients.keys()))
        self._model.update()

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
        if not isinstance(coefficients, dict):
            coefficients = dict(coefficients)
        self._add_constraint(name, coefficients)
        self._model.update()

    def add_columns(self, variables):
        """
        Introduces new variables.

        Parameters
        ----------
        variables: dict or iterable
            dict of dict or an iterable with pairs of variable name and a dict
            of coefficients.
        """
        if not isinstance(variables, dict):
            variables = dict(variables)
        self._model.update()

    def add_rows(self, variables):
        """
        Introduces new variables.

        Parameters
        ----------
        variables: dict or iterable
            dict of dict or an iterable with pairs of variable name and a dict
            of coefficients.
        """
        if not isinstance(variables, dict):
            variables = dict(variables)
        self._model.update()

    def delete_column(self, name):
        """
        Removes a variable from the model.

        Parameters
        ----------
        name: str
            Name of the column variable to be removed.
        """
        self._model.remove(self._variables[name])

    def delete_row(self, name):
        """
        Removes a variable from the model.

        Parameters
        ----------
        name: str
            Name of the row variable to be removed.
        """
        self._model.remove(self._constraints[name])

    def delete_columns(self, variables):
        """
        Removes variables from the model.

        Parameters
        ----------
        variables: iterable
            An iterable that contains the column names as strings.
        """
        map(self.delete_column, variables)

    def delete_rows(self, variables):
        """
        Removes variables from the model.

        Parameters
        ----------
        variables: iterable
            An iterable that contains the row names as strings.
        """
        map(self.delete_row, variables)

    def get_column_names(self):
        """
        Returns
        -------
        list:
            Strings identifying column names.
        """
        return self._variables.keys()

    def get_row_names(self):
        """
        Returns
        -------
        list:
            Strings identifying row names.
        """
        return self._constraints.keys()

    def modify_column_bounds(self, name, constraints):
        """
        Modifies the lower and upper bounds of a particular column.

        Parameters
        ----------
        name: str
            Name of the column variable constraints to be modified.
        constraints: dict or iterable
            Key-value pairs of row names and a pair of lower and upper
            bound.
        """


    def modify_row_bounds(self, name, constraints):
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
        grb_column = self._model.getCol(self._rows[name])
        if not isinstance(constraints, dict):
            constraints = dict(constraints)

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


