from sympy import solve, symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import (print_expression, units, convert_to, Quantity)
from symplyphysics.laws.gravity import free_fall_acceleration_from_height as acceleration

print(f"Formula is:\n{print_expression(acceleration.law)}")

height_above_ground = symbols("height_above_ground")
gravity_constant_value = convert_to(Quantity(units.gravitational_constant),
    units.newton * units.meter**2 / units.kilogram**2).evalf(5)
EARTH_MASS = 5.9722e24  # kilogram
EARTH_RADIUS = 6.371e6  # meter

solved = solve(acceleration.law, acceleration.free_fall_acceleration,
    dict=True)[0][acceleration.free_fall_acceleration]
result_acceleration = solved.subs({
    acceleration.planet_mass: EARTH_MASS,
    acceleration.planet_radius: EARTH_RADIUS,
    acceleration.height_above_surface: height_above_ground,
    acceleration.units.gravitational_constant: gravity_constant_value
})

print(
    f"\nFree fall accelleration function on Earth surface is:\n{print_expression(result_acceleration)}"
)

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
