"""
Braking path via speed and friction force
=========================================

Let an arbitrary object move along the surface at an arbitrary speed. The friction force acts on the
object from the surface. Then the *braking path* will depend on the mass of the object, its speed
and friction force.

**Links:**

#. `Braking distance <https://en.wikipedia.org/wiki/Braking_distance>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as energy_law
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_distance as work_law

braking_path = symbols.distance
"""
The braking path of the object. See :symbols:`distance`.
"""

speed = symbols.speed
"""
The :symbols:`speed` of the object.
"""

friction_force = clone_as_symbol(symbols.force, display_symbol="F_fr", display_latex="F_\\text{fr}")
"""
The friction :symbols:`force` exerted on the object.
"""

mass = symbols.mass
"""
The :symbols:`mass` of the object.
"""

law = Eq(braking_path, mass * speed**2 / (2 * friction_force))
"""
:laws:symbol::

:laws:latex::
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
