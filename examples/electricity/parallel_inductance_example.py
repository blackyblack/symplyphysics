#!/usr/bin/env python3

from sympy import solve, simplify
from symplyphysics import (print_expression, units, Symbol)
from symplyphysics.core.symbols.symbols import tuple_of_symbols
from symplyphysics.laws.electricity import coil_impedance_from_inductivity_and_frequency as coil_impedance_law
from symplyphysics.definitions import admittance_is_inversed_impedance as admittance_def
from symplyphysics.laws.electricity.circuits import admittance_of_parallel_dipoles as parallel_admittance_law

# This example shows how resulting inductivity of 2 parallel coils might be calculated from inductivities of single coils.

# Conditions
## Coils are magnetically decoupled.

inductivity_1 = Symbol("inductivity_1", units.inductance)
inductivity_2 = Symbol("inductivity_2", units.inductance)

# Parallel connection of dipoles summarizes their admittances.
# First find impedances and then admittances
impedance_law = coil_impedance_law.law.subs(coil_impedance_law.coil_impedance,
    admittance_def.dipole_impedance)
admittance_solved = solve([impedance_law, admittance_def.definition],
    admittance_def.dipole_impedance,
    admittance_def.dipole_admittance,
    dict=True)[0][admittance_def.dipole_admittance]
admittance_1 = admittance_solved.subs(coil_impedance_law.coil_inductivity, inductivity_1)
admittance_2 = admittance_solved.subs(coil_impedance_law.coil_inductivity, inductivity_2)
admittances = [admittance_1, admittance_2]

# Then apply parallel admittance law
admittance_symbols = tuple_of_symbols("admittance", units.conductance, len(admittances))
admittances_law = parallel_admittance_law.law.subs(parallel_admittance_law.admittances,
    admittance_symbols).doit()
solved = solve(admittances_law, parallel_admittance_law.parallel_admittance,
    dict=True)[0][parallel_admittance_law.parallel_admittance]
for (from_, to_) in zip(admittance_symbols, admittances):
    solved = solved.subs(from_, to_)
result_admittance = solved

# And finally find resulting inductivity back through the impedance
result_impedance = solve(admittance_def.definition, admittance_def.dipole_impedance,
    dict=True)[0][admittance_def.dipole_impedance].subs(admittance_def.dipole_admittance,
    result_admittance)

result_inductivity = solve(coil_impedance_law.law, coil_impedance_law.coil_inductivity,
    dict=True)[0][coil_impedance_law.coil_inductivity].subs(coil_impedance_law.coil_impedance,
    result_impedance)
result_inductivity = simplify(result_inductivity)
print(f"Inductivity of 2 parallel coils is \n{print_expression(result_inductivity)}")
