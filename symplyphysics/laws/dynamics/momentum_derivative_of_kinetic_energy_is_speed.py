"""
Momentum derivative of kinetic energy is speed
==============================================

The general formula for the kinetic energy of an object features its speed and momentum. This way it can be used
not only in the case of variable mass, but also in the relativistic case.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Kinetic_energy#With_vector_calculus>`__.
"""

from sympy import Eq, Derivative, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_function
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.definitions import momentum_is_mass_times_speed as _momentum_def
from symplyphysics.laws.dynamics import kinetic_energy_via_momentum as _kinetic_energy_law

from symplyphysics.core.experimental.solvers import apply

speed = symbols.speed
"""
The :symbols:`speed` of the object.
"""

momentum = clone_as_function(symbols.momentum, [speed])
"""
The :symbols:`momentum` of the object as a function of :attr:`~speed`.
"""

kinetic_energy = clone_as_function(symbols.kinetic_energy, [momentum(speed)])
"""
The :symbols:`kinetic_energy` of the object as a function of :attr:`~momentum`.
"""

law = Eq(
    Derivative(kinetic_energy(momentum(speed)), momentum(speed)),
    speed,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

_kinetic_energy_eqn = _kinetic_energy_law.law.subs({
    _kinetic_energy_law.kinetic_energy: kinetic_energy(momentum(speed)),
    _kinetic_energy_law.momentum: momentum(speed),
})

_kinetic_energy_diff_momentum_eqn = apply(
    _kinetic_energy_eqn,
    lambda x: x.diff(momentum(speed)),
)

_momentum_eqn = _momentum_def.definition.subs({
    _momentum_def.momentum: momentum(speed),
    _momentum_def.speed: speed,
})

_rhs = solve(
    (_kinetic_energy_diff_momentum_eqn, _momentum_eqn),
    (law.lhs, momentum(speed)),
    dict=True,
)[0][law.lhs]

assert expr_equals(_rhs, law.rhs)


@validate_input(
    kinetic_energy_change_=kinetic_energy,
    momentum_change_=momentum,
)
@validate_output(speed)
def calculate_speed(
    kinetic_energy_change_: Quantity,
    momentum_change_: Quantity,
) -> Quantity:
    kinetic_energy_ = (kinetic_energy_change_ / momentum_change_) * momentum(speed)
    result = law.lhs.subs(
        kinetic_energy(momentum(speed)),
        kinetic_energy_,
    ).doit()
    return Quantity(result)
