#!/usr/bin/env python3

from sympy import pi
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import units, convert_to, quantities
from symplyphysics.laws.quantum_mechanics.harmonic_oscillator import wave_eigenfunctions

values = {
    wave_eigenfunctions.oscillator_mass: convert_to(units.planck_mass, units.kilogram).evalf(),
    wave_eigenfunctions.angular_frequency: 1e10,  # rad/s
}

mode_numbers_ = 0, 1, 2, 3

base_plot = plot(
    title="Quantum harmonic oscillator eigenfunctions",
    xlabel=r"position, $\text{m}$",
    ylabel=r"wave function, $\text{m}^{-1/2}$",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
)

law = wave_eigenfunctions.law.rhs.subs({
    pi: pi.evalf(),
    quantities.hbar: convert_to(units.hbar, units.joule * units.second).evalf(),
})

for mode_number_ in mode_numbers_:
    expr = law.subs(values).subs(wave_eigenfunctions.mode_number, mode_number_)

    sub_plot = plot(
        expr,
        (wave_eigenfunctions.position, -3e-18, 3e-18),  # m
        label=f"$n = {mode_number_}$",
        show=False,
    )
    base_plot.extend(sub_plot)

base_plot.show()
