from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import (
    units, symbols, convert_to, solve, pretty
)
from symplyphysics.laws.gravity import free_fall_acceleration_from_height as acceleration

print("Formula is:\n{}".format(acceleration.print()))

height_above_ground = symbols('height_above_ground')
gravity_constant = convert_to(units.gravitational_constant, units.newton * units.meter**2 / units.kilogram**2).subs(units.newton * units.meter**2 / units.kilogram**2, 1).evalf(5)
earth_mass = 5.976e+24  # kilogram
earth_radius = 6.371e+6 # meter

solved = solve(acceleration.law, acceleration.free_fall_acceleration, dict=True)[0][acceleration.free_fall_acceleration]
result_acceleration = solved.subs({
    acceleration.planet_mass: earth_mass,
    acceleration.planet_radius: earth_radius,
    acceleration.height_above_surface: height_above_ground,
    acceleration.units.gravitational_constant: gravity_constant})

print("\nFree fall accelleration function on Earth surface  is:\n{}".format(
    pretty(result_acceleration, use_unicode=False)))

p1 = plot(
    result_acceleration,
    (height_above_ground, 0, 10000),
    ylim=(9.74, 9.85),
    axis_center=(0.0, 9.75),
    line_color='red',
    title='Acceleration free fall (height)',
    xlabel='height,m',
    ylabel='acceleration,m/s^2',
    label='Acceleration(height)',
    legend=True,
    backend=MatplotlibBackend,
    show=False)
p1.show()
