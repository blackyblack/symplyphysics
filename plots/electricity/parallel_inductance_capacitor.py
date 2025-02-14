#!/usr/bin/env python3

# Description
## Parallel LC circuit might be an object for resonance effect.
## Total LC impedance depends on frequency of supplied energy.
## Impedance has a single peak at the frequency corresponding to Thompson's frequency of this LC.
## This is a resonant frequency of this circuit.
## Also we have here a real LC with some resistive leakage. Resonant frequency stays the same as ideal LC, but resonant peak is not so sharp. The sharpness of resonant peak shows quality factor of LC circuit.

from sympy import Abs, Idx, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import global_index, symbols
from symplyphysics.laws.electricity.circuits import admittance_in_parallel_connection as parallel_admittance_law
from symplyphysics.definitions import admittance_is_inverse_impedance as admittance_def
from symplyphysics.laws.electricity.circuits import capacitor_impedance_from_capacitance_and_frequency as capacitor_impedance
from symplyphysics.laws.electricity.circuits import coil_impedance_via_inductance_and_frequency as coil_impedance
from symplyphysics.laws.electricity.circuits import oscillation_period_of_inductor_capacitor_network as thomsons_formula
from symplyphysics.definitions import period_from_angular_frequency as period_definition

frequency_arg = symbols.angular_frequency

EXAMPLE_INDUCTANCE = 0.00001
EXAMPLE_CAPACITANCE = 0.00001
EXAMPLE_RESISTANCE = 5000000

# Parallel connection of dipoles summarizes their admittances.
# First find impedances and then admittances
L_impedance = coil_impedance.law.rhs.subs({coil_impedance.angular_frequency: frequency_arg})
L_admittance = admittance_def.definition.rhs.subs({
    admittance_def.impedance: L_impedance
}).subs({coil_impedance.inductance: EXAMPLE_INDUCTANCE})

C_impedance = capacitor_impedance.law.rhs.subs(
    {capacitor_impedance.angular_frequency: frequency_arg})
C_admittance = admittance_def.definition.rhs.subs({
    admittance_def.impedance: C_impedance
}).subs({capacitor_impedance.capacitance: EXAMPLE_CAPACITANCE})

R_admittance = admittance_def.definition.rhs.subs({admittance_def.impedance: EXAMPLE_RESISTANCE})

ideal_elements = (L_admittance, C_admittance)
real_elements = (L_admittance, C_admittance, R_admittance)

index_local = Idx("index_local", (1, len(ideal_elements)))
admittance_ideal = parallel_admittance_law.law.rhs.subs(global_index, index_local).doit()
for i, v in enumerate(ideal_elements):
    admittance_ideal = admittance_ideal.subs(parallel_admittance_law.admittance[i + 1], v)

index_local = Idx("index_local", (1, len(real_elements)))
admittance_real = parallel_admittance_law.law.rhs.subs(global_index, index_local).doit()
for i, v in enumerate(real_elements):
    admittance_real = admittance_real.subs(parallel_admittance_law.admittance[i + 1], v)

impedance_from_admittance_law = solve(admittance_def.definition,
    admittance_def.impedance,
    dict=True)[0][admittance_def.impedance]

impedance_ideal = impedance_from_admittance_law.subs({admittance_def.admittance: admittance_ideal})
impedance_ideal_to_plot = Abs(impedance_ideal)

impedance_real = impedance_from_admittance_law.subs({admittance_def.admittance: admittance_real})
impedance_real_to_plot = Abs(impedance_real)

thomsons_period = thomsons_formula.law.rhs.subs({
    thomsons_formula.inductance: EXAMPLE_INDUCTANCE,
    thomsons_formula.capacitance: EXAMPLE_CAPACITANCE
})
frequency_from_period = solve(period_definition.law, period_definition.angular_frequency,
    dict=True)[0][period_definition.angular_frequency]
thomsons_frequency = frequency_from_period.subs({period_definition.period: thomsons_period})

PLOT = plot(impedance_ideal_to_plot,
    (frequency_arg, thomsons_frequency - 1.2, thomsons_frequency + 1.2),
    line_color="blue",
    title="LC impedance",
    label="LC impedance",
    legend=True,
    backend=MatplotlibBackend,
    show=False)

real_impedance = plot(impedance_real_to_plot,
    (frequency_arg, thomsons_frequency - 1.2, thomsons_frequency + 1.2),
    label="RLC impedance",
    line_color="green",
    show=False)
PLOT.append(real_impedance[0])

freq_line = plot((frequency_arg - thomsons_frequency) * 1000000000,
    (frequency_arg, thomsons_frequency - 0.001, thomsons_frequency + 0.001),
    label="frequency = 1 / sqrt(LC)",
    line_color="red",
    show=False)
PLOT.append(freq_line[0])

PLOT.show()
