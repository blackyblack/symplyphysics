from collections import namedtuple
from pytest import fixture, raises
from sympy import diff, sin
from sympy.physics.units import electric_constant
from symplyphysics import (Vector, errors, units, Quantity, assert_equal, prefixes, QuantityVector)
from symplyphysics.core.fields.operators import curl_operator
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.laws.electricity.maxwell_equations import curl_of_magnetic_field_is_conductivity_current_density_and_electric_induction_derivative as current_density_law

## The vector field of magnetic intensity is given. For its distribution, the amplitude is known, which is equal to 1 kiloampere per meter.
## The vector field of electric induction is given. For its distribution, the electric intensity is known, which is equal to 10 kilovolt per meter.
## A point in space in a three-dimensional coordinate system is also given: the x coordinate is 0.493 meter, the y coordinate is 0.5 meter,
## the z coordinate is 0.01249 meter. The time is equal to 1 second. Then components of the conductivity current density vector at a given point
## will be equal to: 877.5825 [ampere / meter^2], 999.922 [ampere / meter^2], 880.917 [ampere / meter^2].

Args = namedtuple("Args", [
    "magnetic_intensity", "electric_induction", "conductivity_current_density", "x", "y", "z",
    "point", "time"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    electric_intensity = Quantity(10 * prefixes.kilo * units.volt / units.meter)
    magnetic_amplitude = Quantity(1 * prefixes.kilo * units.ampere / units.meter)
    time_symbol = current_density_law.time
    time = Quantity(1 * units.second)
    x = Quantity(0.493 * units.meter)
    y = Quantity(0.5 * units.meter)
    z = Quantity(0.01249 * units.meter)
    magnetic_intensity = VectorField(lambda point: [
        sin(point.z / Quantity(1 * units.meter)) * magnetic_amplitude,
        sin(point.x / Quantity(1 * units.meter)) * magnetic_amplitude,
        sin(point.y / Quantity(1 * units.meter)) * magnetic_amplitude
    ])
    electric_induction = VectorField(lambda point: [
        sin(point.x / Quantity(1 * units.meter)) * electric_intensity * electric_constant * sin(
        time_symbol / Quantity(1 * units.second)),
        sin(point.y / Quantity(1 * units.meter)) * electric_intensity * electric_constant * sin(
        time_symbol / Quantity(1 * units.second)),
        sin(point.z / Quantity(1 * units.meter)) * electric_intensity * electric_constant * sin(
        time_symbol / Quantity(1 * units.second))
    ])
    # We need that high precision because electric_constant is very small and resulting electric
    # induction derivative is also very small
    conductivity_current_density = [
        Quantity(877.582561867 * units.ampere / units.meter**2),
        Quantity(999.922000941 * units.ampere / units.meter**2),
        Quantity(880.91701256794 * units.ampere / units.meter**2)
    ]
    return Args(magnetic_intensity=magnetic_intensity,
        electric_induction=electric_induction,
        conductivity_current_density=conductivity_current_density,
        x=x,
        y=y,
        z=z,
        point=(x, y, z),
        time=time)


def test_basic_conductivity_current_density(test_args: Args) -> None:
    result = current_density_law.calculate_conductivity_current_density_at_point(
        test_args.magnetic_intensity, test_args.electric_induction, test_args.point, test_args.time)
    assert_equal(result.components[0], test_args.conductivity_current_density[0])
    assert_equal(result.components[1], test_args.conductivity_current_density[1])
    assert_equal(result.components[2], test_args.conductivity_current_density[2])


def test_basic_magnetic_intensity_rotor(test_args: Args) -> None:
    result = current_density_law.magnetic_intensity_rotor_law(
        test_args.electric_induction,
        Vector(test_args.conductivity_current_density,
        test_args.electric_induction.coordinate_system))
    result_at_point = result.apply(test_args.point)
    result_at_point_quantity = QuantityVector.from_base_vector(result_at_point,
        subs={current_density_law.time: test_args.time})
    rotor_magnetic_intensity = curl_operator(test_args.magnetic_intensity)
    rotor_magnetic_intensity_at_point = rotor_magnetic_intensity.apply(test_args.point)
    assert_equal(result_at_point_quantity.components[0],
        rotor_magnetic_intensity_at_point.components[0])
    assert_equal(result_at_point_quantity.components[1],
        rotor_magnetic_intensity_at_point.components[1])
    assert_equal(result_at_point_quantity.components[2],
        rotor_magnetic_intensity_at_point.components[2])


def test_basic_electric_induction_time_derivative(test_args: Args) -> None:
    result = current_density_law.electric_induction_time_derivative_law(
        test_args.magnetic_intensity,
        Vector(test_args.conductivity_current_density,
        test_args.magnetic_intensity.coordinate_system))
    result_at_point = result.apply(test_args.point)
    electric_induction_at_point = test_args.electric_induction.apply(test_args.point)
    electric_induction_time_derivative = Vector(
        [diff(c, current_density_law.time) for c in electric_induction_at_point.components],
        test_args.electric_induction.coordinate_system)
    electric_induction_at_point_quantity = QuantityVector.from_base_vector(
        electric_induction_time_derivative, subs={current_density_law.time: test_args.time})
    # Increase tolerance due to very small comparables
    assert_equal(result_at_point.components[0],
        electric_induction_at_point_quantity.components[0],
        relative_tolerance=0.1)
    assert_equal(result_at_point.components[1],
        electric_induction_at_point_quantity.components[1],
        relative_tolerance=0.1)
    assert_equal(result_at_point.components[2],
        electric_induction_at_point_quantity.components[2],
        relative_tolerance=0.1)


def test_bad_coordinates(test_args: Args) -> None:
    bad_coordinate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction,
            (bad_coordinate, test_args.y, test_args.z), test_args.time)
    with raises(errors.UnitsError):
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction,
            (test_args.x, bad_coordinate, test_args.z), test_args.time)
    with raises(errors.UnitsError):
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction,
            (test_args.x, test_args.y, bad_coordinate), test_args.time)
    with raises(TypeError):
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction,
            (100, test_args.y, test_args.z), test_args.time)
    with raises(TypeError):
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction,
            (test_args.x, 100, test_args.z), test_args.time)
    with raises(TypeError):
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction,
            (test_args.x, test_args.y, 100), test_args.time)


def test_bad_time(test_args: Args) -> None:
    bad_time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction,
            (test_args.x, test_args.y, test_args.z), bad_time)
    with raises(TypeError):
        current_density_law.calculate_conductivity_current_density_at_point(
            test_args.magnetic_intensity, test_args.electric_induction,
            (test_args.x, test_args.y, test_args.z), 100)
