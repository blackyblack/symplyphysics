from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, units, errors
from symplyphysics.laws.electricity.vector import electric_flux_of_uniform_electric_field as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector

Args = namedtuple("Args", "e a")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e_unit = units.volt / units.meter
    e = QuantityCoordinateVector([3 * e_unit, -3 * e_unit, 0], CARTESIAN)

    a_unit = units.meter**2
    a = QuantityCoordinateVector([0.5 * a_unit, 0.1 * a_unit, 1.2 * a_unit], CARTESIAN)

    return Args(e=e, a=a)


def test_law(test_args: Args) -> None:
    result = law.calculate_electric_flux(test_args.e, test_args.a)
    assert_equal(result, 1.2 * units.volt * units.meter)


def test_bad_electric_field(test_args: Args) -> None:
    eb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_electric_flux(eb_vector, test_args.a)

    eb_scalar = units.volt / units.meter
    with raises(ValueError):
        law.calculate_electric_flux(eb_scalar, test_args.a)

    with raises(TypeError):
        law.calculate_electric_flux(100, test_args.a)
    with raises(TypeError):
        law.calculate_electric_flux([100], test_args.a)


def test_bad_area(test_args: Args) -> None:
    ab_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_electric_flux(test_args.e, ab_vector)

    ab_scalar = units.meter**2
    with raises(ValueError):
        law.calculate_electric_flux(test_args.e, ab_scalar)

    with raises(TypeError):
        law.calculate_electric_flux(test_args.e, 100)
    with raises(TypeError):
        law.calculate_electric_flux(test_args.e, [100])
