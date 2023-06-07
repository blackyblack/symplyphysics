from sympy import solve, symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import (print_expression, units, convert_to)
from symplyphysics.laws.gravity import free_fall_acceleration_from_height as acceleration

print("Formula is:\n{}".format(acceleration.print()))

height_above_ground = symbols("height_above_ground")
gravity_constant_value = convert_to(units.gravitational_constant,
    units.newton * units.meter**2 / units.kilogram**2).evalf(5)
gravity_constant = gravity_constant_value.subs(units.newton * units.meter**2 / units.kilogram**2,
    1).evalf(5)
earth_mass = 5.9722e24  # kilogram
earth_radius = 6.371e6  # meter

solved = solve(acceleration.law, acceleration.free_fall_acceleration,
    dict=True)[0][acceleration.free_fall_acceleration]
result_acceleration = solved.subs({
    acceleration.planet_mass: earth_mass,
    acceleration.planet_radius: earth_radius,
    acceleration.height_above_surface: height_above_ground,
    acceleration.units.gravitational_constant: gravity_constant
})

print("\nFree fall accelleration function on Earth surface is:\n{}".format(
    print_expression(result_acceleration)))

p1 = plot(result_acceleration, (height_above_ground, 0, 10000),
    ylim=(9.74, 9.85),
    axis_center=(0.0, 9.75),
    line_color="red",
    title="Free fall acceleration",
    xlabel="height, m",
    ylabel="acceleration, m/s**2",
    label="Acceleration(height)",
    legend=True,
    backend=MatplotlibBackend,
    show=False)
p1.show()
