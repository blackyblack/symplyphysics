from sympy import solve, pretty, symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import (
    units, convert_to
)
from symplyphysics.laws.gravity import free_fall_acceleration_from_height as acceleration
print("Formula is:\n{}".format(acceleration.print()))
gravity_constant=convert_to(acceleration.gravity_constant, units.newton * units.meter**2 / units.kilogram**2).subs(units.newton * units.meter**2 / units.kilogram**2, 1).evalf(5)
solved = solve(acceleration.law, acceleration.acceleration_free_fall, dict=True)[0][acceleration.acceleration_free_fall]
result_acceleration = solved.subs({acceleration.gravity_constant: gravity_constant,
                                   acceleration.earth_mass: 5.976e+24,
                                   acceleration.earth_radius: 6.371e+6})
print("\nFree fall accelleration function on Earth surface  is:\n{}".format(
    pretty(result_acceleration, use_unicode=False)))
p1 = plot(
    result_acceleration,
    (acceleration.height, 0, 10000),
    ylim=(9.74, 9.85),
    axis_center=(0.0, 9.75),
    line_color='red',
    title='Acceleration free fall (height)',
    xlabel = 'height,m',
    ylabel = 'acceleration,m/s^2',
    label = 'Acceleration(height)',
    legend=True,
    backend = MatplotlibBackend,
    show=False)
p1.show()