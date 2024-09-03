#!/usr/bin/env python3

from sympy import Idx, solve, simplify
from symplyphysics import (print_expression, units, Symbol, global_index)
from symplyphysics.laws.electricity.circuits import coil_impedance_from_inductivity_and_frequency as coil_impedance_law
from symplyphysics.definitions import admittance_is_inverse_impedance as admittance_def
from symplyphysics.laws.electricity.circuits import admittance_in_parallel_connection as parallel_admittance_law

# This example shows how resulting inductivity of 2 parallel coils might be calculated from inductivities of single coils.

# Conditions
## Coils are magnetically decoupled.

inductivity_1 = Symbol("inductivity_1", units.inductance)
inductivity_2 = Symbol("inductivity_2", units.inductance)

# Parallel connection of dipoles summarizes their admittances.
# First find impedances and then admittances
impedance_law = coil_impedance_law.law.subs(coil_impedance_law.coil_impedance,
    admittance_def.impedance)
admittance_solved = solve([impedance_law, admittance_def.definition],
    admittance_def.impedance,
    admittance_def.admittance,
    dict=True)[0][admittance_def.admittance]
admittance_1 = admittance_solved.subs(coil_impedance_law.coil_inductivity, inductivity_1)
admittance_2 = admittance_solved.subs(coil_impedance_law.coil_inductivity, inductivity_2)
admittances = [admittance_1, admittance_2]

# Then apply parallel admittance law
index_local = Idx("index_local", (1, len(admittances)))
admittances_law = parallel_admittance_law.law.subs(global_index, index_local).doit()
for i, v in enumerate(admittances):
    admittances_law = admittances_law.subs(parallel_admittance_law.admittance[i + 1], v)

result_admittance = solve(admittances_law, parallel_admittance_law.total,
    dict=True)[0][parallel_admittance_law.total]

# And finally find resulting inductivity back through the impedance
result_impedance = solve(admittance_def.definition, admittance_def.impedance,
    dict=True)[0][admittance_def.impedance].subs(admittance_def.admittance, result_admittance)

result_inductivity = solve(coil_impedance_law.law, coil_impedance_law.coil_inductivity,
    dict=True)[0][coil_impedance_law.coil_inductivity].subs(coil_impedance_law.coil_impedance,
    result_impedance)
result_inductivity = simplify(result_inductivity)
print(f"Inductivity of 2 parallel coils is \n{print_expression(result_inductivity)}")
