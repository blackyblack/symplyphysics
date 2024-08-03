from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, QuantityVector)
from symplyphysics.core.fields.operators import curl_operator
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.vectors.arithmetics import diff_cartesian_vector, scale_vector
from symplyphysics.laws.electricity.maxwell_equations import derivative_of_magnetic_induction_in_time_is_rotor_of_electric_intensity as faradays_law

Args = namedtuple("Args",
    ["electric_intensity", "magnetic_induction", "x", "y", "z", "point", "time"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    electric_intensity_unit = Quantity(units.volt / units.meter)
    magnetic_induction_amplitude = Quantity(10 * units.tesla)
    time = Quantity(1 * units.second)
    x = Quantity(0.493 * units.meter)
    y = Quantity(0.5 * units.meter)
    z = Quantity(0.01249 * units.meter)
    magnetic_induction = VectorField(
        [0, 0, faradays_law.time / Quantity(1 * units.second) * magnetic_induction_amplitude])
    electric_intensity = VectorField(lambda point: [(magnetic_induction_amplitude * point.y /
        2) / units.meter / units.tesla * electric_intensity_unit,
        (-1 * magnetic_induction_amplitude * point.x / 2) / units.meter / units.tesla *
        electric_intensity_unit,
        Quantity(0)])
    return Args(electric_intensity=electric_intensity,
        magnetic_induction=magnetic_induction,
        x=x,
        y=y,
        z=z,
        point=(x, y, z),
        time=time)


def test_basic_magnetic_induction_derivative(test_args: Args) -> None:
    result = faradays_law.calculate_magnetic_induction_derivative_at_point(
        test_args.electric_intensity, test_args.point, test_args.time)
    magnetic_induction_space = test_args.magnetic_induction.apply_to_basis()
    magnetic_induction_derivative = diff_cartesian_vector(magnetic_induction_space,
        faradays_law.time)
    magnetic_induction_derivative = scale_vector(-1, magnetic_induction_derivative)
    assert_equal(result.components[0], magnetic_induction_derivative.components[0])
    assert_equal(result.components[1], magnetic_induction_derivative.components[1])
    assert_equal(result.components[2], magnetic_induction_derivative.components[2])


def test_basic_electric_intensity_curl(test_args: Args) -> None:
    result = faradays_law.electric_intensity_curl_law(test_args.magnetic_induction)
    result_at_point = result.apply(test_args.point)
    result_at_point_quantity = QuantityVector.from_base_vector(result_at_point,
        subs={faradays_law.time: test_args.time})
    rotor_electric_intensity = curl_operator(test_args.electric_intensity)
    rotor_electric_intensity_at_point = rotor_electric_intensity.apply(test_args.point)
    assert_equal(result_at_point_quantity.components[0],
        rotor_electric_intensity_at_point.components[0])
    assert_equal(result_at_point_quantity.components[1],
        rotor_electric_intensity_at_point.components[1])
    assert_equal(result_at_point_quantity.components[2],
        rotor_electric_intensity_at_point.components[2])


def test_bad_coordinates(test_args: Args) -> None:
    bad_coordinate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            (bad_coordinate, test_args.y, test_args.z), test_args.time)
    with raises(errors.UnitsError):
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            (test_args.x, bad_coordinate, test_args.z), test_args.time)
    with raises(errors.UnitsError):
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            (test_args.x, test_args.y, bad_coordinate), test_args.time)
    with raises(TypeError):
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            (100, test_args.y, test_args.z), test_args.time)
    with raises(TypeError):
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            (test_args.x, 100, test_args.z), test_args.time)
    with raises(TypeError):
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            (test_args.x, test_args.y, 100), test_args.time)


def test_bad_time(test_args: Args) -> None:
    bad_time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            (test_args.x, test_args.y, test_args.z), bad_time)
    with raises(TypeError):
        faradays_law.calculate_magnetic_induction_derivative_at_point(test_args.electric_intensity,
            (test_args.x, test_args.y, test_args.z), 100)
