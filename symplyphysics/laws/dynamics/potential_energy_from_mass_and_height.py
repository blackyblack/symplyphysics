"""
Potential energy from mass and height
=====================================

The potential energy is a form of energy a body possesses because of its position relative to
other objects.

**Conditions:**

#. Acceleration due to gravity is constant throughout the region of the body's movement.

**Links:**

#. `Wikipedia, second formula <https://en.wikipedia.org/wiki/Potential_energy#Overview>`__.

..
    TODO: add note that it is specifically gravitational potential energy => move to `gravity`?
"""

from sympy import Eq, solve, symbols as sym_symbols
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities,
    clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, CoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector

from symplyphysics.definitions import (
    mechanical_energy_is_kinetic_and_potential_energy as _energy_def,)
from symplyphysics.laws.conservation import (
    initial_mechanical_energy_equals_final_mechanical_energy as _const_energy_law)
from symplyphysics.laws.dynamics import total_work_is_change_in_kinetic_energy as _work_energy_law
from symplyphysics.laws.dynamics.vector import (
    acceleration_from_force_vector as _newtons_second_law,
    mechanical_work_from_force_and_displacement as _work_law,
)

potential_energy = symbols.potential_energy
"""
The :symbols:`potential_energy` of the body.
"""

mass = symbols.mass
"""
The :symbols:`mass` of the body.
"""

height = symbols.height
"""
The elevation from ground level. See :symbols:`height`
"""

law = Eq(potential_energy, mass * quantities.acceleration_due_to_gravity * height)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

# 1. Find the work done by the body to move in the (Earth's) gravitational field.

# Note that the potential energy the body possesses is equal (by magnitude) to the amount of work
# it can do against the gravity force, so we can use the expression for work.

_free_fall_acceleration = CoordinateVector(
    [0, 0, -1 * quantities.acceleration_due_to_gravity],
    CARTESIAN,
)

_gravity_force = solve_for_vector(
    _newtons_second_law.law,
    _newtons_second_law.force,
).subs({
    _newtons_second_law.mass: mass,
    _newtons_second_law.acceleration: _free_fall_acceleration,
})

_dx, _dy = sym_symbols("dx:y", real=True)

# Displacement from initial point A to final point B:
_displacement = CoordinateVector([_dx, _dy, height], CARTESIAN)

# The work that the object would perform against the gravitational field to move from A to B.
# Note that in order to use this law, we presume that the gravity force (and therefore the
# acceleration due to gravity) is constant throughout the region of the object's movement.
_work_expr = _work_law.law.rhs.subs({
    _work_law.force: _gravity_force,
    _work_law.displacement: _displacement
})

# 2. Show that `U = -W`

_k0 = _work_energy_law.kinetic_energy(_work_energy_law.time_before)
_k1 = _work_energy_law.kinetic_energy(_work_energy_law.time_after)

_u0 = clone_as_symbol(potential_energy, subscript="0")
_u1 = clone_as_symbol(potential_energy, subscript="1")

_e0 = _energy_def.definition.rhs.subs({
    _energy_def.kinetic_energy: _k0,
    _energy_def.potential_energy: _u0,
})

_e1 = _energy_def.definition.rhs.subs({
    _energy_def.kinetic_energy: _k1,
    _energy_def.potential_energy: _u1,
})

_total_energy_eqn = _const_energy_law.law.subs({
    _const_energy_law.mechanical_energy(_const_energy_law.initial_time): _e0,
    _const_energy_law.mechanical_energy(_const_energy_law.final_time): _e1,
})

_work_energy_eqn = _work_energy_law.law.subs({
    _work_energy_law.work: _work_expr,
})

# The "potential energy" used in the law is actually defined as the the difference in potential
# energy (or, more precisely, the gravitational potential) after and before the movement.
_potential_energy_eqn = Eq(
    potential_energy,
    _u1 - _u0,
)

_potential_energy_derived = solve(
    (_total_energy_eqn, _work_energy_eqn, _potential_energy_eqn),
    (potential_energy, _k0, _u0),
    dict=True,
)[0][potential_energy]

assert expr_equals(_potential_energy_derived, law.rhs)


@validate_input(body_mass_=mass, height_=height)
@validate_output(potential_energy)
def calculate_potential_energy(body_mass_: Quantity, height_: Quantity) -> Quantity:
    result_energy_expr = solve(law, potential_energy, dict=True)[0][potential_energy]
    result_expr = result_energy_expr.subs({mass: body_mass_, height: height_})
    return Quantity(result_expr)
