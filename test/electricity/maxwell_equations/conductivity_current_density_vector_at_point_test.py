from collections import namedtuple
from pytest import fixture, raises
from sympy import sin
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
    prefixes
)
from symplyphysics.core.fields.vector_field import VectorField
from sympy.physics.units import electric_constant

from symplyphysics.laws.electricity.maxwell_equations import conductivity_current_density_vector_at_point as current_density_law

## The vector of magnetic intensity is given. For its distribution, the amplitude is known, which is equal to 1 kiloampere per meter.
## The vector of electrical induction is given. For its distribution, the electric intensity is known, which is equal to 10 kilovolt per meter.
## A point in space in a three-dimensional coordinate system is also given: the x coordinate is 0.493 meter, the y coordinate is 0.5 meter,
## the z coordinate is 0.01249 meter. The time is equal to 1 second. Then components of the conductivity current density vector at a given point
## will be equal to: 877 [ampere / meter^2], 999 [ampere / meter^2], 881 [ampere / meter^2].

Args = namedtuple("Args", ["magnetic_intensity", "electrical_induction", "x", "y", "z", "time"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    electric_intensity = Quantity(10 * prefixes.kilo * units.volt / units.meter)
    mag_amplitude = Quantity(1 * prefixes.kilo * units.ampere / units.meter)
    time = Quantity(1 * units.second)
    x = Quantity(0.493 * units.meter)
    y = Quantity(0.5 * units.meter)
    z = Quantity(0.01249 * units.meter)
    magnetic_intensity = VectorField(lambda point: [sin(point.z / Quantity(1 * units.meter)) * mag_amplitude,
                                                    sin(point.x / Quantity(1 * units.meter)) * mag_amplitude,
                                                    sin(point.y / Quantity(1 * units.meter)) * mag_amplitude])
    electrical_induction = VectorField(lambda point: [sin(point.x / Quantity(1 * units.meter)) * electric_intensity * electric_constant * sin(time / Quantity(1 * units.second)),
                                                      sin(point.y / Quantity(1 * units.meter)) * electric_intensity * electric_constant * sin(time / Quantity(1 * units.second)),
                                                      sin(point.z / Quantity(1 * units.meter)) * electric_intensity * electric_constant * sin(time / Quantity(1 * units.second))])
    return Args(magnetic_intensity=magnetic_intensity,
        electrical_induction=electrical_induction,
        x=x,
        y=y,
        z=z,
        time=time)


def test_basic_conductivity_current_density(test_args: Args) -> None:
    result = current_density_law.calculate_conductivity_current_density_at_point(test_args.magnetic_intensity, test_args.electrical_induction, (test_args.x, test_args.y, test_args.z), test_args.time)
    assert_equal(result.components[0], 877 * units.ampere / units.meter**2)
    assert_equal(result.components[1], 999 * units.ampere / units.meter**2)
    assert_equal(result.components[2], 881 * units.ampere / units.meter**2)


def test_bad_coordinates(test_args: Args) -> None:
    bad_coordinate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_density_law.calculate_conductivity_current_density_at_point(test_args.magnetic_intensity, test_args.electrical_induction, (bad_coordinate, test_args.y, test_args.z), test_args.time)
    with raises(errors.UnitsError):
        current_density_law.calculate_conductivity_current_density_at_point(test_args.magnetic_intensity, test_args.electrical_induction, (test_args.x, bad_coordinate, test_args.z), test_args.time)
    with raises(errors.UnitsError):
        current_density_law.calculate_conductivity_current_density_at_point(test_args.magnetic_intensity, test_args.electrical_induction, (test_args.x, test_args.y, bad_coordinate), test_args.time)
    with raises(TypeError):
        current_density_law.calculate_conductivity_current_density_at_point(test_args.magnetic_intensity, test_args.electrical_induction, (100, test_args.y, test_args.z), test_args.time)
    with raises(TypeError):
        current_density_law.calculate_conductivity_current_density_at_point(test_args.magnetic_intensity, test_args.electrical_induction, (test_args.x, 100, test_args.z), test_args.time)
    with raises(TypeError):
        current_density_law.calculate_conductivity_current_density_at_point(test_args.magnetic_intensity, test_args.electrical_induction, (test_args.x, test_args.y, 100), test_args.time)


def test_bad_time(test_args: Args) -> None:
    bad_time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_density_law.calculate_conductivity_current_density_at_point(test_args.magnetic_intensity, test_args.electrical_induction, (test_args.x, test_args.y, test_args.z), bad_time)
    with raises(TypeError):
        current_density_law.calculate_conductivity_current_density_at_point(test_args.magnetic_intensity, test_args.electrical_induction, (test_args.x, test_args.y, test_args.z), 100)
