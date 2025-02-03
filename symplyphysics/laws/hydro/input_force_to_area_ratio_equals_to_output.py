"""
Force to area ratio in hydraulic press
======================================

If both vertically positioned cylinders of communicating vessels are closed with
pistons, then with the help of external forces applied to the pistons, a large pressure
can be created in the liquid, many times exceeding the hydrostatic pressure at any point
in the system. If the pistons have different areas, then different forces act on them
from the liquid side. The same modulus, but oppositely directed external forces must be
applied to the pistons to keep the system in balance.


**Conditions:**

#. This ratio is performed only in an ideal hydraulic press, i.e. one in which there is
   no friction.

**Links:**

#. `Physics LibreTexts, formula 53.1.2 <https://phys.libretexts.org/Courses/Prince_Georges_Community_College/General_Physics_I%3A_Classical_Mechanics/53%3A_Hydraulics_and_Pneumatics/53.01%3A_Hydraulics-_The_Hydraulic_Press>`__.

..
    TODO: rename file
"""

from sympy import Eq, solve, dsolve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.laws.dynamics import pressure_from_force_and_area as pressure_law
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.hydro import inner_pressure_is_constant as constant_pressure_law

input_force = clone_as_symbol(symbols.force, subscript="1")
"""
:symbols:`force` acting on the first piston.
"""

input_area = clone_as_symbol(symbols.area, subscript="1")
"""
:symbols:`area` of the first piston.
"""

output_force = clone_as_symbol(symbols.force, subscript="2")
"""
:symbols:`force` acting on the second piston.
"""

output_area = clone_as_symbol(symbols.area, subscript="2")
"""
:symbols:`area` of the second piston.
"""

law = Eq(input_force / input_area, output_force / output_area)
"""
:laws:symbol::

:laws:latex::
"""

# TODO prefix variables used in proof with underscore

_pressure_input = pressure_law.law.rhs.subs({
    pressure_law.force: input_force,
    pressure_law.area: input_area
})

_pressure_output = pressure_law.law.rhs.subs({
    pressure_law.force: output_force,
    pressure_law.area: output_area,
})

## If the pistons are in equilibrium, then the pressures _pressure_input and _pressure_output are equal
_dsolved = dsolve(constant_pressure_law.law,
    constant_pressure_law.inner_pressure(constant_pressure_law.time))
_dsolved_input = _dsolved.subs(constant_pressure_law.inner_pressure(constant_pressure_law.time),
    _pressure_input)
_dsolved_output = _dsolved.subs(constant_pressure_law.inner_pressure(constant_pressure_law.time),
    _pressure_output)
_solved_input = solve([_dsolved_input, _dsolved_output], (_pressure_input, "C1"),
    dict=True)[0][_pressure_input]
_pressure_equation = Eq(_pressure_input, _solved_input)

assert expr_equals(law.rhs, _pressure_equation.rhs)
assert expr_equals(law.lhs, _pressure_equation.lhs)


@validate_input(input_force_=input_force,
    input_area_=input_area,
    output_forces_area_=output_area)
@validate_output(output_force)
def calculate_output_force(input_force_: Quantity, input_area_: Quantity,
    output_forces_area_: Quantity) -> Quantity:
    result_expr = solve(law, output_force, dict=True)[0][output_force]
    result_force = result_expr.subs({
        input_force: input_force_,
        input_area: input_area_,
        output_area: output_forces_area_,
    })
    return Quantity(result_force)
