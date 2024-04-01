#!/usr/bin/env python3

from collections import namedtuple
from sympy import pi
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.laws.waves import displacement_in_interfering_waves as interference_law

# Description
## Plot the interference of two waves of unit amplitude and angular wavenumber at zero time
## depending on the position for different values of the phase shift between the waves.

wave = interference_law.law.rhs.subs({
    interference_law.amplitude: 1,
    interference_law.angular_wavenumber: 1,
    interference_law.time: 0,
})

print(f"Expression of the interfering waves:\n{print_expression(wave)}\n")

Data = namedtuple("Data", "phi label")

plot_data_ = (
    Data(phi=0, label="0"),
    Data(phi=pi / 3, label=r"\frac{\pi}{3}"),
    Data(phi=pi / 2, label=r"\frac{\pi}{2}"),
    Data(phi=2 * pi / 3, label=r"\frac{2\pi}{3}"),
    Data(phi=pi, label=r"\pi"),
)

p = plot(
    title="Interference of waves",
    xlabel="position, m",
    ylabel="displacement, m",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
    annotations=None,
)

for plot_datum in plot_data_:
    wave_subs = wave.subs(interference_law.phase_shift, plot_datum.phi)
    sub_p = plot(
        wave_subs,
        (interference_law.position, -5, 5),
        label=f"$\\phi = {plot_datum.label}$",
        show=False,
    )
    p.append(sub_p[0])

p.show()
