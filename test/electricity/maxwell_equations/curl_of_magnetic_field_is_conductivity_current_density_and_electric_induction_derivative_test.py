from collections import namedtuple
from pytest import fixture, raises
from sympy import sin
from symplyphysics import errors, units, Quantity, units, quantities
from symplyphysics.laws.electricity.maxwell_equations import (
    curl_of_magnetic_field_is_conductivity_current_density_and_electric_induction_derivative as
    current_density_law)

from symplyphysics.core.experimental.coordinate_systems import (CARTESIAN, CoordinateVector,
    AppliedPoint, QuantityCoordinateVector)
from symplyphysics.core.experimental.approx import assert_equal_vectors

## The vector field of magnetic intensity is given. For its distribution, the amplitude is known, which is equal to 1 kiloampere per meter.
## The vector field of electric induction is given. For its distribution, the electric intensity is known, which is equal to 10 kilovolt per meter.
## A point in space in a three-dimensional coordinate system is also given: the x coordinate is 0.493 meter, the y coordinate is 0.5 meter,
## the z coordinate is 0.01249 meter. The time is equal to 1 second. Then components of the conductivity current density vector at a given point
## will be equal to: 877.5825 [ampere / meter^2], 999.922 [ampere / meter^2], 880.917 [ampere / meter^2].

Args = namedtuple("Args", [
    "magnetic_intensity",
    "electric_induction",
    "conductivity_current_density",
    "point",
    "time",
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    electric_intensity = Quantity(10 * units.kilo * units.volt / units.meter)
    magnetic_amplitude = Quantity(1 * units.kilo * units.ampere / units.meter)

    x, y, z = CARTESIAN.base_scalars

    magnetic_intensity = CoordinateVector([
        sin(z / units.meter) * magnetic_amplitude,
        sin(x / units.meter) * magnetic_amplitude,
        sin(y / units.meter) * magnetic_amplitude
    ], CARTESIAN)

    t = current_density_law.time

    electric_constant = quantities.vacuum_permittivity
    electric_induction = CoordinateVector([
        sin(x / units.meter) * electric_intensity * electric_constant * sin(t / units.second),
        sin(y / units.meter) * electric_intensity * electric_constant * sin(t / units.second),
        sin(z / units.meter) * electric_intensity * electric_constant * sin(t / units.second)
    ], CARTESIAN)

    # We need that high precision because electric_constant is very small and resulting electric
    # induction derivative is also very small
    conductivity_current_density = QuantityCoordinateVector([
        Quantity(877.582561867 * units.ampere / units.meter**2),
        Quantity(999.922000941 * units.ampere / units.meter**2),
        Quantity(880.91701256794 * units.ampere / units.meter**2)
    ], CARTESIAN)

    point = AppliedPoint([
        Quantity(0.493 * units.meter),
        Quantity(0.5 * units.meter),
        Quantity(0.01249 * units.meter),
    ], CARTESIAN)

    time = units.second

    return Args(
        magnetic_intensity=magnetic_intensity,
        electric_induction=electric_induction,
        conductivity_current_density=conductivity_current_density,
        point=point,
        time=time,
    )


def test_basic_conductivity_current_density(test_args: Args) -> None:
    result = current_density_law.calculate_conductivity_current_density_at_point(
        test_args.magnetic_intensity, test_args.electric_induction, test_args.point, test_args.time)
    assert_equal_vectors(result, test_args.conductivity_current_density)


def test_bad_coordinates(test_args: Args) -> None:
    bad_coordinate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bad_point = AppliedPoint([bad_coordinate, 0, 0], CARTESIAN)
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction, bad_point, test_args.time)
    with raises(TypeError):
        bad_point = AppliedPoint([100, 0, 0], CARTESIAN)
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction, bad_point, test_args.time)

    with raises(errors.UnitsError):
        bad_point = AppliedPoint([0, bad_coordinate, 0], CARTESIAN)
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction, bad_point, test_args.time)
    with raises(TypeError):
        bad_point = AppliedPoint([0, 100, 0], CARTESIAN)
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction, bad_point, test_args.time)

    with raises(errors.UnitsError):
        bad_point = AppliedPoint([0, 0, bad_coordinate], CARTESIAN)
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction, bad_point, test_args.time)
    with raises(TypeError):
        bad_point = AppliedPoint([0, 0, 100], CARTESIAN)
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction, bad_point, test_args.time)


def test_bad_time(test_args: Args) -> None:
    bad_time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction, test_args.point, bad_time)
    with raises(TypeError):
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction, test_args.point, 100)
