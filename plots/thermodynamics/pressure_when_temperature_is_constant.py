#!/usr/bin/env python3
from sympy import solve, symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.laws.thermodynamics import pressure_and_volume_in_isothermal_process as boyles_law

print(f"Formula is:\n{print_expression(boyles_law.law)}")

volume = symbols("volume")

solved = solve(
    boyles_law.law,
    boyles_law.final_pressure,
    dict=True,
)[0][boyles_law.final_pressure]

result_pressure = solved.subs({
    boyles_law.initial_pressure: 1,
    boyles_law.initial_volume: 1,
    boyles_law.final_volume: volume
})

print(f"Pressure function is:\n{print_expression(result_pressure)}")

p1 = plot(result_pressure, (volume, 0.01, 1),
    title="Pressure(Volume)",
    xlabel="Volume",
    ylabel="Pressure",
    backend=MatplotlibBackend,
    show=False)

p1.show()
