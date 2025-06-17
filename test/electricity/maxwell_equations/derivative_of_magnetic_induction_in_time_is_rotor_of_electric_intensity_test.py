from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.electricity.maxwell_equations import (
    derivative_of_magnetic_induction_in_time_is_rotor_of_electric_intensity as faradays_law)

from symplyphysics.core.experimental.vectors import vector_diff
from symplyphysics.core.experimental.coordinate_systems import (CARTESIAN, CoordinateVector,
    AppliedPoint, QuantityCoordinateVector)
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", ["electric_intensity", "magnetic_induction", "point", "time"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    electric_intensity_unit = Quantity(units.volt / units.meter)

    magnetic_induction_amplitude = Quantity(10 * units.tesla)

    time = Quantity(1 * units.second)

    point = AppliedPoint([
        Quantity(0.493 * units.meter),
        Quantity(0.5 * units.meter),
        Quantity(0.01249 * units.meter),
    ], CARTESIAN)

    magnetic_induction = CoordinateVector(
        [0, 0, 1 * faradays_law.time / Quantity(1 * units.second) * magnetic_induction_amplitude],
        CARTESIAN)

    x, y, _ = CARTESIAN.base_scalars

    magnetic_factor = ((magnetic_induction_amplitude * electric_intensity_unit) /
        (2 * units.meter * units.tesla))

    electric_intensity = CoordinateVector(
        [magnetic_factor * y, -1 * magnetic_factor * x, 0],
        CARTESIAN,
    )

    return Args(
        electric_intensity=electric_intensity,
        magnetic_induction=magnetic_induction,
        point=point,
        time=time,
    )


def test_basic_magnetic_induction_derivative(test_args: Args) -> None:
    result = faradays_law.calculate_magnetic_induction_derivative_at_point(
        test_args.electric_intensity, test_args.point, test_args.time)

    expected = vector_diff(test_args.magnetic_induction, faradays_law.time)
    expected = CoordinateVector.from_expr(expected)
    expected = QuantityCoordinateVector(expected.components, expected.system, expected.point)

    assert_equal_vectors(result, expected)


def test_bad_coordinates(test_args: Args) -> None:
    bad_coordinate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bad_point = AppliedPoint((bad_coordinate, 0, 0), CARTESIAN)
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            bad_point, test_args.time)
    with raises(TypeError):
        bad_point = AppliedPoint((100, 0, 0), CARTESIAN)
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            bad_point, test_args.time)

    with raises(errors.UnitsError):
        bad_point = AppliedPoint((0, bad_coordinate, 0), CARTESIAN)
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            bad_point, test_args.time)
    with raises(TypeError):
        bad_point = AppliedPoint((0, 100, 0), CARTESIAN)
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            bad_point, test_args.time)

    with raises(errors.UnitsError):
        bad_point = AppliedPoint((0, 0, bad_coordinate), CARTESIAN)
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            bad_point, test_args.time)
    with raises(TypeError):
        bad_point = AppliedPoint((0, 0, 100), CARTESIAN)
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            bad_point, test_args.time)


def test_bad_time(test_args: Args) -> None:
    bad_time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            test_args.point, bad_time)
    with raises(TypeError):
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            test_args.point, 100)
