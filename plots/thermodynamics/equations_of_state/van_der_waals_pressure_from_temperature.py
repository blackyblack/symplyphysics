#!/usr/bin/env python3

from sympy import symbols, Eq, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import equation as van_der_waals_law
from symplyphysics.laws.quantities import quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law

volume_values = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]  # liters

temperature = symbols("temperature")
pressure = symbols("pressure")
volume = symbols("volume")
amount_of_substance = symbols("amount_of_substance")

parameter_a = symbols("parameter_a")
parameter_b = symbols("parameter_b")

molar_volume = solve(molar_qty_law.law, molar_qty_law.molar_quantity)[0].subs({
    molar_qty_law.extensive_quantity: volume,
    molar_qty_law.amount_of_substance: amount_of_substance,
})

state_equation = van_der_waals_law.law.subs({
    van_der_waals_law.pressure: pressure,
    van_der_waals_law.temperature: temperature,
    van_der_waals_law.molar_volume: molar_volume,
    van_der_waals_law.attractive_forces_parameter: parameter_a,
    van_der_waals_law.excluded_volume_parameter: parameter_b
})
pressure_value = solve(state_equation, pressure, dict=True)[0][pressure]
answer = Eq(pressure, pressure_value)
print(f"Total equation:\n{print_expression(answer)}")

pressure_to_plots = pressure_value.subs({
    parameter_a: 0.188,  # [Pa * m^6 / mole^2]
    parameter_b: 4.532 * 1E-5,  # [m^3 / mole]
    amount_of_substance: 1,  # [moles]
})
pressure_to_plots = evaluate_expression(pressure_to_plots)

base_plot = plot(title="The Van der Waals state equation",
    xlabel=r"$T, K$",
    ylabel=r"$P, Pa$",
    backend=MatplotlibBackend,
    legend=True,
    show=False)

for volume_value in volume_values:
    # Convert from liters to m^3
    pressure_to_subplot = pressure_to_plots.subs({volume: volume_value * 1E-3})
    subplot = plot(
        pressure_to_subplot,
        (temperature, 100, 500),  # [K]
        label=f"Volume is {volume_value}",
        show=False)
    base_plot.append(subplot[0])

base_plot.show()
