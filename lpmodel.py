#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
====================
LP Solver Interfaces
====================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
:Date:
    2011-03-28
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    lpmodel.py
"""


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
        name: str:
            Identifier for the column.
        coefficients: dict or iterable:
            key-value pairs of row names and their coefficients.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def add_row(self, name, coefficients):
        """
        Introduces a new variable.

        Parameters
        ----------
        name: str:
            Identifier for the row.
        coefficients: dict or iterable:
            key-value pairs of column names and their coefficients.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def add_columns(self, variables):
        """
        Introduces new variables.

        Parameters
        ----------
        variables: dict or iterable:
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
        variables: dict or iterable:
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
        name: str:
            Name of the column variable to be removed.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def delete_row(self, name):
        """
        Removes a variable from the model.

        Parameters
        ----------
        name: str:
            Name of the row variable to be removed.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def delete_columns(self, variables):
        """
        Removes variables from the model.

        Parameters
        ----------
        variables: iterable:
            An iterable that contains the column names as strings.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def delete_rows(self, variables):
        """
        Removes variables from the model.

        Parameters
        ----------
        variables: iterable:
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

    def modify_column_constraints(self, name, constraints):
        """
        Modifies the lower and upper bounds of a particular column.

        Parameters
        ----------
        name: str:
            Name of the column variable constraints to be modified.
        constraints: dict or iterable:
            Key-value pairs of row names and a pair of lower and upper
            bound.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def modify_row_constraints(self, name, constraints):
        """
        Modifies the lower and upper bounds of a particular row.

        Parameters
        ----------
        name: str:
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
        variables: iterable:
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
    """

    def __init__(self, name=""):
        """
        Parameters
        ----------
        name: str:
            Optional name of the model.
        """
        LPModelFacade.__init__(self)
        self._model = gurobipy.Model(name)
        self._columns = dict()
        self._rows = dict()

    def add_column(self, name, coefficients):
        """
        Introduces a new variable.

        Parameters
        ----------
        name: str:
            Identifier for the column.
        coefficients: dict or iterable:
            key-value pairs of row names and their coefficients.
        """
        coeffs = list()
        rows = list()
        if not isinstance(coefficients, dict):
            coefficients = dict(coefficients)
        for (key, value) in coefficients.iteritems():
            rows.append(self._rows[key])
            coeffs.append(value)
        self._model.addConstr(gurobipy.LinExpr(coeffs, rows), gurobipy.EQUAL,
                0.0, name)
        self._columns[name] = set(self._model.getConstrs()).difference(
                set(self._columns.values()))
        self._model.update()

    def add_row(self, name, coefficients):
        """
        Introduces a new variable.

        Parameters
        ----------
        name: str:
            Identifier for the row.
        coefficients: dict or iterable:
            key-value pairs of column names and their coefficients.
        """
        coeffs = list()
        columns = list()
        if not isinstance(coefficients, dict):
            coefficients = dict(coefficients)
        for (key, value) in coefficients.iteritems():
            columns.append(self._rows[key])
            coeffs.append(value)
        column = gurobipy.Column(coeffs, columns)
        self._rows[name] = self._model.addVar(name=name, column=column)
        self._model.update()

    def add_columns(self, variables):
        """
        Introduces new variables.

        Parameters
        ----------
        variables: dict or iterable:
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
        variables: dict or iterable:
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
        name: str:
            Name of the column variable to be removed.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def delete_row(self, name):
        """
        Removes a variable from the model.

        Parameters
        ----------
        name: str:
            Name of the row variable to be removed.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def delete_columns(self, variables):
        """
        Removes variables from the model.

        Parameters
        ----------
        variables: iterable:
            An iterable that contains the column names as strings.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def delete_rows(self, variables):
        """
        Removes variables from the model.

        Parameters
        ----------
        variables: iterable:
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

    def modify_column_constraints(self, name, constraints):
        """
        Modifies the lower and upper bounds of a particular column.

        Parameters
        ----------
        name: str:
            Name of the column variable constraints to be modified.
        constraints: dict or iterable:
            Key-value pairs of row names and a pair of lower and upper
            bound.
        """
        raise NotImplementedError("abstract base class, subclass to expose an "\
                "interface to a specific LP solver")

    def modify_row_constraints(self, name, constraints):
        """
        Modifies the lower and upper bounds of a particular row.

        Parameters
        ----------
        name: str:
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
        variables: iterable:
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


