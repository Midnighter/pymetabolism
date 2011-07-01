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


import errno


class Error(StandardError):
    """
    An error for all exceptions that occur in the usage of the pymetabolism
    package.
    """

    def __init__(self, msg, num=1, *args, **kw_args):
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
        >>> err = pymetabolism.Error("It's too %s outside!", "rainy")
        >>> print(err)
        It's too rainy outside!
        >>> print(err.errno)
        1
        """
        StandardError.__init__(self, *args, **kw_args)
        self.args = (msg,) + args
        self.errno = num
        self.strerror = msg % args

    def __str__(self):
        """
        Returns
        -------
        str:
            Simply returns the formatted string.
        """
        return self.strerror

