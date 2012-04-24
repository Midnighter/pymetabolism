#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
======================================
Network degree distribution generators
======================================

:Author:
    Alexandra Grigore
:Date:
    2012-03-23
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    degree_distribution.py
"""

import pylab as P
import random
import numpy


def powerlaw_sequence(n,expo=2.0):
    """
    Returns a power law degree distribution.
    
    Parameters
    ----------
    n: int
        The length of the power law sequence
    expo: float
        The exponent of the power law sequence
    """
    return numpy.random.zipf(expo,n)
    
def normal_sequence(n,mu,sigma):
    """
    Returns a normal degree distribution.
    
    Parameters
    ----------
    n: int
        The length of the normal degree sequence
    mu: int
        The mean of the normal distribution
    sigma: int
        The standard deviation of the normal distribution
    """
    return numpy.random.normal(mu,sigma,n)
    
#def plot_powerlaw(sequence,bins=10):
#    """
#    Plots the power law degree distribution in a loglog plot
#    """
#    h,be = numpy.histogram(sequence,bins)
#    P.figure()
#    P.loglog(h)
#    P.show()
#    
#def plot_normal(sequence,bins=10):
#    """
#    Plots the normal degree distribution
#    """
#    P.figure()    
#    h, bins, patches = P.hist(sequence, bins, normed=1, histtype='stepfilled')
#    P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
#    P.show()
    

