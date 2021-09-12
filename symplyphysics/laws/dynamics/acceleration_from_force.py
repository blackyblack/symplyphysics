from sympy import symbols, Eq, pretty, solve
from sympy.physics.units import Quantity
import sympy.physics.units as phy_units
from symplyphysics.quantity_decorator import validate_input, validate_output
import symplyphysics.expr_to_quantity

# Description
## Newton's second law: a = F / m

force, mass, acceleration = symbols('force mass acceleration')
law = Eq(acceleration, force / mass)

def print():
    return pretty(law, use_unicode=False)

@validate_input(mass_=phy_units.mass, acceleration_=phy_units.acceleration)
@validate_output(phy_units.force)
def calculate_force(mass_: Quantity, acceleration_: Quantity) -> Quantity:
    result_expr = solve(law.subs(mass, mass_).subs(acceleration, acceleration_))[0]
    return symplyphysics.expr_to_quantity.convert(result_expr, 'force')
