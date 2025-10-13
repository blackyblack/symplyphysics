"""
Acceleration from force and mass (vector)
=========================================

Newton's second law of motion states that in an inertial frame of reference, the acceleration of
a body in motion is proportional to the net force exerted on it, and its mass is the constant of
proportionality.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Newton's_laws_of_motion#Second_law>`__.
"""

from sympy import Eq, Expr
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol

from symplyphysics.core.experimental.solvers import solve_for_vector, vector_equals
from symplyphysics.core.experimental.vectors import (clone_as_vector_symbol,
    clone_as_vector_function)
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

from symplyphysics.definitions.vector import (
    momentum_is_mass_times_velocity_vector as _momentum_def,
    acceleration_is_velocity_derivative as _acceleration_def,
)
from symplyphysics.laws.dynamics.vector import force_is_derivative_of_momentum as _newtons_law

mass = clone_as_symbol(symbols.mass, positive=True)
"""
:symbols:`mass` of the object.
"""

acceleration = clone_as_vector_symbol(symbols.acceleration)
"""
Vector of the :symbols:`acceleration` of the object.
"""

force = clone_as_vector_symbol(symbols.force)
"""
Vector of the :symbols:`force` exerted on the object.
"""

law = Eq(acceleration, force / mass)
"""
:laws:symbol::

:laws:latex::
"""

# Derive this law from law of force and momentum
# Condition: mass is constant

_time = _newtons_law.time
_velocity = clone_as_vector_function(symbols.speed, (_time,))

_momentum_expr = solve_for_vector(_momentum_def.law, _momentum_def.momentum).subs({
    _momentum_def.mass: mass,
    _momentum_def.velocity: _velocity(_time),
})

_newtons_law_subs = _newtons_law.law.subs({
    _newtons_law.force(_time): force,
    _newtons_law.momentum(_time): _momentum_expr,
})

_force_derived = solve_for_vector(_newtons_law_subs, force)

_acceleration_def_subs = _acceleration_def.law.subs({
    _acceleration_def.acceleration(_acceleration_def.time): acceleration,
    _acceleration_def.velocity(_acceleration_def.time): _velocity(_time),
})

_force_expected = solve_for_vector(law, force).subs(
    _acceleration_def_subs.lhs,
    _acceleration_def_subs.rhs,
)

assert vector_equals(_force_derived, _force_expected)


@validate_input(mass_=mass, acceleration_=acceleration)
@validate_output(force)
def calculate_force(
    mass_: Quantity,
    acceleration_: QuantityCoordinateVector,
) -> Expr:
    force_expr = solve_for_vector(law, force)
    force_value = force_expr.subs({
        mass: mass_,
        acceleration: acceleration_,
    })

    return QuantityCoordinateVector.from_expr(force_value)


@validate_input(mass_=mass, force_=force)
@validate_output(acceleration)
def calculate_acceleration(
    mass_: Quantity,
    force_: QuantityCoordinateVector,
) -> Expr:
    acceleration_expr = solve_for_vector(law, acceleration)
    acceleration_value = acceleration_expr.subs({
        mass: mass_,
        force: force_,
    })

    return QuantityCoordinateVector.from_expr(acceleration_value)
