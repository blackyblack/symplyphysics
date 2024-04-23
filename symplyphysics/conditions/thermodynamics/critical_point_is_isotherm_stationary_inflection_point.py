from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Symbol,
    Function,
)

# Description
## Critical point (in the thermodynamic sense) is such values of volume, pressure, and temperature at
## which only one phase exists and at the vicinity of which the physical properties of the phases of the
## substance change dramatically. Algebraically, the critical point is the stationary inflection point
## of the isothermal pressure-volume dependency line.

# Conditions: (dp/dV)_T = 0, (d**2(p)/dV**2)_T = 0
## p - pressure
## V - volume
## T - temperature
## (d/dV)_T - partial derivative with respect to volume at constant temperature

# Note
## - These equations need to be solved together with the equation of state.

pressure = Function("pressure", units.pressure)
volume = Symbol("volume", units.volume)

inflection_point_condition = Eq(Derivative(pressure(volume), volume, 1), 0)
flat_tangent_condition = Eq(Derivative(pressure(volume), volume, 2), 0)

conditions = [inflection_point_condition, flat_tangent_condition]
