#!/usr/bin/env python3

from sympy import symbols, solve, Rational
from sympy.plotting import plot, plot_parametric
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import convert_to, units, quantities
from symplyphysics.laws.thermodynamics import (
    pressure_and_volume_in_adiabatic_process as adiabatic_law,)
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation as ideal_gas_law

ideal_gas_law_subs = ideal_gas_law.law.subs({
    quantities.molar_gas_constant: units.molar_gas_constant,
})

pressure, volume, temperature = symbols("pressure volume temperature", positive=True)

ideal_gas_eqn = ideal_gas_law_subs.subs({
    ideal_gas_law.amount_of_substance:
        1,
    units.molar_gas_constant:
    convert_to(
    units.molar_gas_constant,
    units.joule / (units.kelvin * units.mol),
    ),
    ideal_gas_law.pressure:
        pressure,
    ideal_gas_law.volume:
        volume,
    ideal_gas_law.temperature:
        temperature,
})

pressure_expr = solve(ideal_gas_eqn, pressure)[0]

base_plot = plot(
    title="Fundamental thermodynamic processes of an ideal monatomic gas",
    xlabel=r"volume, $\text{m}^3$",
    ylabel=r"pressure, $\text{Pa}$",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
)

# Plot guiding isotherms

for temperature_, label in zip((400, 800, 1200), "lower middle upper".split()):
    pressure_expr_ = pressure_expr.subs(temperature, temperature_)
    isotherm_plot = plot(
        pressure_expr_,
        (volume, 1, 5),
        # ideally the label should go to the right end of the corresponding line
        # although I reckon this is only possible to do via matplotlib
        label=f"$T = {temperature_} \\, \\text{{K}}$ ({label})",
        line_color="yellow",
        show=False,
    )
    base_plot.extend(isotherm_plot)

# Let the prossesses start at V = 1.5 m**3 and T = 800 K

INITIAL_VOLUME = 1.5  # m**3

INITIAL_TEMPERATURE = 800  # K
FINAL_TEMPERATURE = 400  # K

initial_pressure = pressure_expr.subs({
    volume: INITIAL_VOLUME,
    temperature: INITIAL_TEMPERATURE,
})

# Plot adiabate

# For monatomic gas with 3 translational degrees of freedom
adiabatic_index = Rational(5, 3)

adiabate_eqn = adiabatic_law.adiabatic_condition.subs({
    adiabatic_law.adiabatic_index: adiabatic_index,
    adiabatic_law.initial_pressure: initial_pressure,
    adiabatic_law.initial_volume: INITIAL_VOLUME,
    adiabatic_law.final_pressure: pressure,
    adiabatic_law.final_volume: volume,
})

final_volume = solve(
    [
    adiabate_eqn,
    ideal_gas_eqn.subs(temperature, FINAL_TEMPERATURE),
    ],
    (pressure, volume),
    dict=True,
)[0][volume]

adiabate_pressure_expr = solve(adiabate_eqn, pressure)[0]

adiabate_plot = plot(
    adiabate_pressure_expr,
    (volume, INITIAL_VOLUME, final_volume),
    label=r"adiabate, $Q = 0$",
    line_color="blue",
    show=False,
)
base_plot.extend(adiabate_plot)

# Plot isobar

isobar_plot = plot(
    initial_pressure,
    (volume, INITIAL_VOLUME, final_volume),
    label=r"isobar, $p = \text{const}$",
    line_color="green",
    show=False,
)
base_plot.extend(isobar_plot)

# Plot isotherm

isotherm_plot = plot(
    pressure_expr.subs(temperature, INITIAL_TEMPERATURE),
    (volume, INITIAL_VOLUME, final_volume),
    label=r"isotherm, $T = \text{const}$",
    line_color="red",
    show=False,
)
base_plot.extend(isotherm_plot)

# Plot isochore

isochore_final_pressure = pressure_expr.subs({
    volume: INITIAL_VOLUME,
    temperature: FINAL_TEMPERATURE,
})

isochore_plot = plot_parametric(
    (INITIAL_VOLUME, pressure),
    (pressure, initial_pressure, isochore_final_pressure),
    label=r"isochore, $V = \text{const}$",
    line_color="violet",
    show=False,
)
base_plot.extend(isochore_plot)

base_plot.show()
