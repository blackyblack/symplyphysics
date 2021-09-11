from sympy import symbols, Eq, pretty, solve
import sympy.physics.units as phy_units
from sympy.physics.units import Quantity
from sympy.physics.units.systems.si import SI

# Description
## Newton's second law: a = F / m

force, mass, acceleration = symbols('force mass acceleration')
law = Eq(acceleration, force / mass)

def print():
    return pretty(law, use_unicode=False)

def calculate_force(mass_: Quantity, acceleration_: Quantity) -> Quantity:
    dimsys_SI = SI.get_dimension_system()
    assert dimsys_SI.equivalent_dims(mass_.dimension, phy_units.mass)
    assert dimsys_SI.equivalent_dims(acceleration_.dimension, phy_units.acceleration)

    result_expr = solve(law.subs(mass, mass_).subs(acceleration, acceleration_))[0]
    quantity_scale = SI._collect_factor_and_dimension(result_expr)
    result_force = Quantity('force')
    dimsys_SI.set_quantity_dimension(result_force, quantity_scale[1])
    dimsys_SI.set_quantity_scale_factor(result_force, quantity_scale[0])
    assert dimsys_SI.equivalent_dims(result_force.dimension, phy_units.force)
    return result_force