#!/usr/bin/env python3

from sympy import  symbols, plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, units
from symplyphysics.laws.electricity import electric_field_due_to_point_charge as electric_field

# Description
## Reference to example (net_field_is_zero_from_two_point_charges.py)[../../../examples/electricity/net_field_is_zero_from_two_point_charges.py]
## Plot the net field of two point charges along the x-axis.

COULOMB_CONSTANT = 8.99e9  # N * m^2 * C^(-2)

CHARGE_1 = 2.1e-8  # C
CHARGE_2 = -4.0 * CHARGE_1

POSITION_OF_CHARGE_1 = 0.2  # m
POSITION_OF_CHARGE_2 = 0.7  # m

position = symbols("position")

field_of_charge_1 = electric_field.law.rhs.subs({
    electric_field.point_charge: CHARGE_1,
    electric_field.distance: position - POSITION_OF_CHARGE_1,
    units.coulomb_constant: COULOMB_CONSTANT,
})
field_of_charge_2 = electric_field.law.rhs.subs({
    electric_field.point_charge: CHARGE_2,
    electric_field.distance: position - POSITION_OF_CHARGE_2,
    units.coulomb_constant: COULOMB_CONSTANT,
})
net_field = field_of_charge_1 + field_of_charge_2

print(f"Net field expression:\n{print_expression(net_field)}")

plot_of_total_field = plot(
    net_field,
    (position, -0.5, 1),
    title="Net field of two point charges",
    label="total field",
    line_color="blue",
    xlabel="x, m",
    ylabel="E, V/m",
    ylim=(-3e6, 3e6),
    legend=True,
    backend=MatplotlibBackend,
    show=False,
)

plot_of_total_field.show()
