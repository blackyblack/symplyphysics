#!/usr/bin/env python3

from sympy import solve
from symplyphysics import (units, convert_to, Quantity, prefixes, Symbol)
from symplyphysics.laws.electricity import coil_impedance_from_inductivity_and_frequency as coil_impedance_law
from symplyphysics.definitions import admittance_is_inversed_impedance as admittance_def
from symplyphysics.laws.electricity.circuits import admittance_of_parallel_dipoles as parallel_admittance_law

# This example shows how resulting inductivity of 2 parallel coils might be calculated from inductivities of single coils.

# Conditions
## Coils are magnetically decoupled.

inductivity_1 = Symbol("inductivity_1", units.inductance)
inductivity_2 = Symbol("inductivity_2", units.inductance)
w = Symbol("frequency", units.frequency)

# Parallel connection of dipoles summarizes their admittances.
# First find impedances and then admittances
impedance_1 = solve(coil_impedance_law.law, coil_impedance_law.coil_impedance, dict=True)[0][coil_impedance_law.coil_impedance].subs({
    coil_impedance_law.coil_inductivity: inductivity_1,
    coil_impedance_law.circular_frequency: w
})
impedance_2 = solve(coil_impedance_law.law, coil_impedance_law.coil_impedance, dict=True)[0][coil_impedance_law.coil_impedance].subs({
    coil_impedance_law.coil_inductivity: inductivity_2,
    coil_impedance_law.circular_frequency: w
})

admittance_1 = solve(admittance_def.definition, admittance_def.dipole_admittance, dict=True)[0][admittance_def.dipole_admittance].subs({
    admittance_def.dipole_impedance: impedance_1
})
admittance_2 = solve(admittance_def.definition, admittance_def.dipole_admittance, dict=True)[0][admittance_def.dipole_admittance].subs({
    admittance_def.dipole_impedance: impedance_2
})

admittances = {admittance_1, admittance_2}

# Then apply parallel admittance law
result_admittance = solve(parallel_admittance_law.law, parallel_admittance_law.parallel_admittance, dict=True)[0][parallel_admittance_law.parallel_admittance].subs({
    parallel_admittance_law.admittances: admittances
})

# And finally find resulting inductivity back through the impedance
result_impedance = solve(admittance_def.definition, admittance_def.dipole_impedance, dict=True)[0][admittance_def.dipole_impedance].subs({
    admittance_def.dipole_admittance: result_admittance
})

result_inductivity = "1"

print(f"Inductivity of 2 parallel coils is {result_inductivity}")
