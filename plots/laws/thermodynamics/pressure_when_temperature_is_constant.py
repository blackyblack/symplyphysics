#!/usr/bin/env python3
from sympy import solve, symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.laws.thermodynamics import temperature_is_constant as boyles_law

print("Formula is:\n{}".format(boyles_law.print()))
volume = symbols("volume")
solved = solve(boyles_law.law, boyles_law.pressure_end, dict=True)[0][boyles_law.pressure_end]
result_pressure = solved.subs({
    boyles_law.pressure_start: 1,
    boyles_law.volume_start: 1,
    boyles_law.volume_end: volume
})

print("Pressure function is:\n{}".format(print_expression(result_pressure)))

p1 = plot(result_pressure, (volume, 0.01, 1),
    title="Pressure(Volume)",
    xlabel="Volume",
    ylabel="Pressure",
    backend=MatplotlibBackend,
    show=False)

p1.show()
