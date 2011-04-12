#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
====================
Metabolic Components
====================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
    Nikolaus Sonnenschein
:Date:
    2011-04-07
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    metabolism.py
"""


class BasicMetabolicComponent(object):
    """
    Base class for all components of a metabolic model system.

    Notes
    -----
    Each component is uniquely identified and stored by its class and name.
    Instantiating the same class with the same name will simply yield the
    original instance. Metabolic components can only be compared for equality or
    inequality.

    Examples
    --------

    >>> comp = BasicMetabolicComponent("M1")
    >>> comp_2 = BasicMetabolicComponent("M1")
    >>> comp == comp_2
    True

    """

    _counter = 1
    _memory = dict()

    def __new__(cls, name="", *args, **kw_args):
        """
        Ensures the unique instance policy of all metabolic components.
        """
        if cls._memory.has_key((cls, name)):
            return cls._memory[(cls, name)]
        else:
            return object.__new__(cls)

    def __init__(self, name=""):
        """
        Parameters
        ----------
        name: str (optional)
            A string uniquely identifying one component among its class.
        """
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        object.__init__(self)
        self._index = self.__class__._counter
        self.__class__._counter += 1
        if name:
            self.name = name
        else:
            self.name = "%s_%d" % (self.__class__.__name__, self._index)
        self.__class__._memory[(self.__class__, self.name)] = self

    def __str__(self):
        return self.name

    def __repr__(self):
        return "<%s.%s, %d>" % (self.__module__, self.__class__.__name__, self._index)

    def __lt__(self, other):
        raise NotImplementedError

    def __le__(self, other):
        raise NotImplementedError

    def __eq__(self, other):
        if not isinstance(other, BasicMetabolicComponent):
            raise TypeError("invalid comparison between %s and %s" %
                    (repr(self), repr(other)))
        elif self.__class__ == other.__class__:
            return self._index == other._index
        else:
            return False

    def __ne__(self, other):
        if not isinstance(other, BasicMetabolicComponent):
            raise TypeError("invalid comparison between %s and %s" %
                    (repr(self), repr(other)))
        elif self.__class__ == other.__class__:
            return self._index != other._index
        else:
            return True

    def __gt__(self, other):
        raise NotImplementedError

    def __ge__(self, other):
        raise NotImplementedError

    def __cmp__(self, other):
        raise NotImplementedError


class SimpleCompound(BasicMetabolicComponent):
    """
    The simplest form of representing a metabolic compound - just its class and
    name.
    """

    def __init__(self, name=""):
        """
        Parameters
        ----------
        name: str (optional)
            A string uniquely identifying the compound among its class.
        """
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        BasicMetabolicComponent.__init__(self, name=name)


class SimpleReaction(BasicMetabolicComponent):
    """
    The simplest form of representing a biochemical reaction - just its class
    and name.

    Notes
    -----
    This is useful in metabolic networks where substrates and products are
    determined by the topology. Stoichiometric information is then stored on the
    links themselves. The only other information stored about the reaction is
    its reversibility.
    """

    def __init__(self, name="", reversible=False):
        """
        Parameters
        ----------
        name: str (optional)
            A string uniquely identifying the reaction among its class.
        reversible: bool (optional)
            Reversibility information of the reaction.
        """
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        BasicMetabolicComponent.__init__(self, name=name)
        self.reversible = reversible

class SBMLCompartment(BasicMetabolicComponent):
    """
    A cellular compartment as defined per SBML standard.
    """

    def __init__(self, name, outside=None, constant=True, suffix="",
            spatial_dimensions=None, size=None, units=None):
        """
        Parameters
        ----------
        name: str
            A string uniquely identifying the compartment among its class.
        outside: str
            The name of the compartment that surrounds this one.
        constant: bool (optional)
            Determines whether the size attribute is allowed to change during
            model simulation.
        suffix: str (optional)
            A string appended to compounds for input/output.
        spatial_dimensions: int (optional)
            From 0 to 3, normal models have three dimensions.
        size: float (optional)
            The magnitude of the spatial_dimension in units.
        units: str (optional)
            A string identifying the unit in which size is measured.

        Notes
        -----
        The constant attribute is so far only kept for compatibility with SBML,
        it's not actually required. This behaviour may change in future.
        """
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        BasicMetabolicComponent.__init__(self, name=name)
        if outside:
            self.outside = SBMLCompartment(outside)
        else:
            self.outside = None
        self.constant = bool(constant)
        self.suffix = str(suffix)
        self.spatial_dimensions = spatial_dimensions
        self.size = size
        self.units = units


class SBMLCompound(BasicMetabolicComponent):
    """
    A molecular compound as defined per SBML standard.
    """

    def __init__(self, identifier, formula=None, in_chl=None, in_chl_key=None,
            smiles=None, charge=None, mass=None):
        """
        Parameters
        ----------
        identifier: str
            A string uniquely identifying the compound among its class.
        formula: str (optional)
            Molecular formula as a simple string, e.g., C6H12O6.
        in_chl: str (optional)
            An IUPAC compliant identifier in InChl format.
        in_chl_key: int (optional)
            A hashed key of the InChl string.
        smiles: str (optional)
            A SMILES representation of the compound.
        charge: int (optional)
            Electric charge on the compound (may be pH dependent).
        mass: float (optional)
            A unit-less magnitude determining the mass of the compound.
        """
        if self.__class__._memory.has_key((self.__class__, identifier)):
            return
        BasicMetabolicComponent.__init__(self, name=identifier)
        self.identifier = identifier
        self.formula = formula
        self.in_chl = in_chl
        self.in_chl_key = in_chl_key
        self.smiles = smiles
        try:
            self.charge = int(charge)
        except (ValueError, TypeError):
            self.charge = None
        try:
            self.mass = float(mass)
        except (ValueError, TypeError):
            self.mass = None

    def __contains__(self, element):
        """
        Checks for the existance of an atomic element in the compound.
        """
        raise NotImplementedError


class SBMLCompartmentCompound(BasicMetabolicComponent):
    """
    A compartment specific compound as defined per SBML standard.

    TODO:
        Make direct attribute reference to self.compound work by defining
        __getattr__ (previous implementations failed on recursion)
    """

    def __new__(cls, compound, compartment):
        """
        Ensures the unique instance policy of all metabolic components.
        """
        if compartment.suffix:
            name = compound.identifier + compartment.suffix
        else:
            name = "%s(%s)" % (compound.identifier, compartment.name)
        if cls._memory.has_key((cls, name)):
            return cls._memory[(cls, name)]
        else:
            return BasicMetabolicComponent.__new__(cls, name, compound,
                    compartment)

    def __init__(self, compound, compartment):
        """
        Parameters
        ----------
        compound: SBMLCompound
            An instance of SBMLCompound that is then attached to a compartment.
        compartment: SBMLCompartment
            An instance of SBMLCompartment in which the compound is located.
        """
        if compartment.suffix:
            name = compound.identifier + compartment.suffix
        else:
            name = "%s(%s)" % (compound.identifier, compartment.name)
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        BasicMetabolicComponent.__init__(self, name=name)
        self.compound = compound
        self.compartment = compartment

    def __getattr__(self, attr):
        return self.compound.__getattribute__(attr)

class SBMLReaction(BasicMetabolicComponent):
    """
    A biochemical reaction as defined per SBML standard.
    """

    def __init__(self, identifier, substrates, products, reversible=False,
            synonyms=None, rate_constant=None, lower_bound=None,
            upper_bound=None, objective_coefficient=None, flux_value=None,
            reduced_cost=None):
        """
        Parameters
        ----------
        identifier: str
            A string uniquely identifying the reaction among its class.
        substrates: dict
            A map from the reaction educts to the aboslute value of their
            stoichiometric factors in the reaction.
        products: dict
            A map from the reaction products to the aboslute value of their
            stoichiometric factors in the reaction.
        reversible: bool (optional)
            Whether this reaction is known to occur in both directions in an
            organism.
        synonyms: str (optional)
            Additional identifiers of the reaction.
        rate_constant: float (optional)
            Unit-less specifier of the rate of the reaction at model conditions.
        """
        if self.__class__._memory.has_key((self.__class__, identifier)):
            return
        BasicMetabolicComponent.__init__(self, name=identifier)
        self.substrates = substrates
        self.products = products
        self.reversible = bool(reversible)
        self.synonyms = synonyms
        try:
            self.rate_constant = float(rate_constant)
        except (ValueError, TypeError):
            self.rate_constant = None
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.objective_coefficient = objective_coefficient
        self.flux_value = flux_value
        self._consistency_check()

    def __contains__(self, compound):
        """
        Parameters
        ----------
        compound: SBMLCompound
            A compound instance whose participation in the reaction is tested.
        """
        return (compound in self.substrates or compound in self.products)

    def __len__(self):
        return len(self.substrates) + len(self.products)

    def __str__(self):
        """
        Returns
        -------
        str:
            A string representation of the reaction, e.g., '2 A + 4 B -> 1 C'
            or '2 A + 4 B <=> 1 C' for a reversible reaction.
        """

        def util(compounds):
            for cmpd in compounds:
                yield str(abs(self.stoichiometric_coeff(cmpd)))
                yield str(cmpd)
                if not (cmpd == compounds[-1]):
                    yield "+"

        rxn = [e for e in util(self.substrates.keys())]
        if self.reversible:
            rxn.append("<=>")
        else:
            rxn.append("->")
        rxn.extend([e for e in util(self.products.keys())])
        return " ".join(rxn)

    def compounds(self):
        """
        Returns
        -------
        iterator:
            An iterator over all compounds partaking in the reaction.
        """
        return (cmpd for cmpd in self.substrates.keys() + self.products.keys())

    def stoichiometric_coeff(self, compound):
        """
        Parameters
        ----------
        compound: SBMLCompound
            A compound instance whose stoichiometric coefficient is sought for.

        Returns
        -------
        float:
            The stoichiometric coefficient of a compound in the reaction.
            Coefficients of substrates are negative.

        Exceptions
        ----------
        KeyError:
            In case compound is not part of the reaction.
        """
        if compound in self.substrates:
            return -self.substrates[compound]
        elif compound in self.products:
            return self.products[compound]
        else:
            raise KeyError("'%s' is not participating in reaction '%s'"\
                % (str(compound), self.identifier))

    def _consistency_check(self):
        """
        Asserts some basic consistency of the SBMLReaction instance.
        With enough meta data (SBMLCompound formula, charge, or mass)
        stoichiometric balancing is checked.

        Exceptions
        ----------
        AssertionError:
            In case any of the given conditions are not true.
        """
        # elemental balancing
        if all(cmpd.formula for cmpd in self.substrates.keys() +
                self.products.keys()):
            pass # not implemented yet
        # mass balancing
        if all(cmpd.mass for cmpd in self.substrates.keys() +
                self.products.keys()):
            assert sum(cmpd.mass * coeff for (cmpd, coeff) in
                    self.substrates.iteritems()) == sum(cmpd.mass * coeff
                    for (cmpd, coeff) in self.products.iteritems()),\
                    "There is a mass imbalance in reaction '%s'" % self.identifier
        # charge balancing
        if all(cmpd.charge for cmpd in self.substrates.keys() +
                self.products.keys()):
            assert sum(cmpd.charge * coeff for (cmpd, coeff) in
                    self.substrates.iteritems()) == sum(cmpd.charge * coeff
                    for (cmpd, coeff) in self.products.iteritems()),\
                    "There is a charge imbalance in reaction '%s'" % self.identifier

    def is_substrate(self, compound):
        """
        Parameters
        ----------
        compound: SBMLCompound
            A compound instance whose status as educt or product is queried.
        """
        return compound in self.substrates

