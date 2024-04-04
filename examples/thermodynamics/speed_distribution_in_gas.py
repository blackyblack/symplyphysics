#!/usr/bin/env python3

from sympy import solve, Symbol, Eq, symbols
from symplyphysics import print_expression, Quantity, prefixes, units, convert_to, quantities
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution

# Description
## In oxygen, whose molar mass is 0.032 kg/mol, at room temperature, what fraction of the molecules
## have speeds in the interval 599 to 601 m/s?

molar_mass = symbols("molar_mass")
