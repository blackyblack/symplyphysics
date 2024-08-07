"""
Braking path via speed and friction force
=========================================

Let an arbitrary object move along the surface at an arbitrary speed. The friction force acts on the
object from the surface. Then the *braking path* will depend on the mass of the object, its speed
and friction force.

..
    TODO Rename file
"""

from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as energy_law
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_move as work_law

braking_path = Symbol("braking_path", units.length)
"""
The braking path of the object.

Symbol:
    :code:`l`
"""

speed = Symbol("speed", units.speed)
"""
The speed of the object.

Symbol:
    :code:`v`
"""

friction_force = clone_symbol(symbols.dynamics.force, "friction_force")
r"""
The friction :attr:`~symplyphysics.symbols.dynamics.force` exerted on the object.

Symbol:
    :code:`F_fr`

Latex:
    :math:`F_\text{fr}`
"""

mass = symbols.basic.mass
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the object.

Symbol:
    :code:`m`
"""

law = Eq(braking_path, mass * speed**2 / (2 * friction_force))
r"""
:code:`l = m * v^2 / (2 * F_fr)`

Latex:
    .. math::
        l = \frac{m v^2}{2 F_\text{fr}}
"""

# This law might be derived via "kinetic_energy_from_mass_and_velocity" law and
# "mechanical_work_from_force_and_move" law.

_energy_law_applied = energy_law.law.subs({energy_law.mass: mass, energy_law.speed: speed})
_energy_derived = solve(_energy_law_applied, energy_law.kinetic_energy,
    dict=True)[0][energy_law.kinetic_energy]

_work_law_applied = work_law.law.subs({
    work_law.force: friction_force,
    work_law.work: _energy_derived
})
_distance_derived = solve(_work_law_applied, work_law.distance, dict=True)[0][work_law.distance]

# Check if derived distance is same as declared.
assert expr_equals(_distance_derived, law.rhs)


@validate_input(mass_=mass, velocity_=speed, friction_force_=friction_force)
@validate_output(braking_path)
def calculate_braking_path(mass_: Quantity, velocity_: Quantity,
    friction_force_: Quantity) -> Quantity:
    result_braking_path_expr = solve(law, braking_path, dict=True)[0][braking_path]
    result_expr = result_braking_path_expr.subs({
        mass: mass_,
        speed: velocity_,
        friction_force: friction_force_
    })
    return Quantity(result_expr)
