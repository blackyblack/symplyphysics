#!/usr/bin/env python3

from sympy import symbols, Eq, solve, simplify, pi
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.laws.optics import refraction_angle_from_environments as refraction_law

# Description
## Plot the dependency of refraction angle on incidence angle for different pairs of environments

AIR_REFRACTIVE_INDEX = 1.0003
WATER_REFRACTIVE_INDEX = 1.333
BENZENE_REFRACTIVE_INDEX = 1.501
DIAMOND_REFRACTIVE_INDEX = 2.417


