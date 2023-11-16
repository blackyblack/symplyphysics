#!/usr/bin/env python3

# Description
## Parallel LC circuit might be an object for resonans effect.
## Total LC impedance depends on frequency of supplied energy.
## Impedance have a single peak at the frequency corresponding to Tompson's frequency of this LC.
## This is a resonant frequency of this circuit.

from sympy import Abs, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import Symbol
from symplyphysics.laws.electricity.circuits import admittance_of_parallel_dipoles as parallel_admittance_law
from symplyphysics.definitions import admittance_is_inversed_impedance as admittance_def
from symplyphysics.laws.electricity import capacitor_impedance_from_capacitance_and_frequency as capacitor_impedance
from symplyphysics.laws.electricity import coil_impedance_from_inductivity_and_frequency as coil_impedance
from symplyphysics.laws.electricity.circuits import oscillation_period_for_capacitor_inductor_node as tompsons_formula
from symplyphysics.laws.kinematic import period_from_angular_frequency as period_definition

frequency_arg = Symbol("frequency_arg")

EXAMPLE_CAPACITANCE = 0.00001
EXAMPLE_INDUCTANCE = 0.00001

# Parallel connection of dipoles summarizes their admittances.
# First find impedances and then admittances
L_impedance = coil_impedance.law.rhs.subs({coil_impedance.circular_frequency : frequency_arg})
L_admittance = admittance_def.definition.rhs.subs({admittance_def.dipole_impedance: L_impedance})

C_impedance = capacitor_impedance.law.rhs.subs({capacitor_impedance.circular_frequency : frequency_arg})
C_admittance = admittance_def.definition.rhs.subs({admittance_def.dipole_impedance: C_impedance})

total_impedance_from_admittance_law = solve(admittance_def.definition, admittance_def.dipole_impedance, dict=True)[0][admittance_def.dipole_impedance]

total_impedance_complex = total_impedance_from_admittance_law.subs({admittance_def.dipole_admittance: (L_admittance + C_admittance)})

total_impedance_abs_to_plot = Abs(total_impedance_complex).subs({coil_impedance.coil_inductivity: EXAMPLE_INDUCTANCE, capacitor_impedance.capacitor_capacitance: EXAMPLE_CAPACITANCE})

tomsons_period = tompsons_formula.law.rhs.subs({tompsons_formula.inductance: EXAMPLE_INDUCTANCE, tompsons_formula.capacitance: EXAMPLE_CAPACITANCE})
frequency_from_period = solve(period_definition.law, period_definition.circular_frequency, dict = True)[0][period_definition.circular_frequency]
tomsons_frequency = frequency_from_period.subs({period_definition.period: tomsons_period})


PLOT = plot(total_impedance_abs_to_plot, (frequency_arg, tomsons_frequency - 1.5, tomsons_frequency + 1.5),
    line_color="blue",
    title="LC impedance",
    label="LC impedance",
    legend=True,
    annotations={},
    backend=MatplotlibBackend,
    show=False)

freq_line = plot((frequency_arg - tomsons_frequency) * 1000000000,  (frequency_arg, tomsons_frequency - 0.001, tomsons_frequency + 0.001),
    label="frequency = 1 / sqrt(LC)",
    line_color="red",    
    show=False)
PLOT.append(freq_line[0])

PLOT.show()
