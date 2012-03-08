#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
==========
Exceptions
==========

:Author:
    Moritz Emanuel Beber
:Date:
    2011-04-13
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    errors.py
"""


class PyMetabolismError(StandardError):
    """
    An error for all exceptions that occur in the usage of the pymetabolism
    package.
    """

    def __init__(self, msg, *args, **kw_args):
        """
        Parameters
        ----------
        msg: str
            An unformatted string, i.e., it may contain multiple string format
            markers.

        Notes
        -----
        A variable number of arguments may be passed. They will all be used to
        format msg. So take care that the number and type of additional
        arguments matches the format markers in msg.

        Examples
        --------
        >>> err = PyMetabolismError("It's too %s outside!", "rainy")
        >>> print(err)
        It's too rainy outside!
        >>> print(err.errno)
        1
        """
        super(PyMetabolismError, self).__init__()
        self.args = (msg,) + args
        self.errno = kw_args.get("errno", 1)
        try:
            self.strerror = msg % kw_args
        except TypeError:
            self.strerror = msg % args

    def __str__(self):
        """
        Returns
        -------
        str:
            Simply returns the formatted string.
        """
        return self.strerror

