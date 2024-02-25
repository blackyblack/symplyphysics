from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.hydro import pressure_of_liquid_in_vessel_moving_horizontally as pressure_law

# Description
## Consider a vessel 1 meter high with water moving horizontally with an acceleration of 1 [meter / second^2].
## Then the water pressure in the vessel will be 9827 pascal.

Args = namedtuple("Args", ["density_liquid", "acceleration", "height"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    density_liquid = Quantity(997 * (units.kilogram / units.meter**3))
    acceleration = Quantity(1 * units.meter / units.second**2)
    height = Quantity(1 * units.meter)

    return Args(density_liquid=density_liquid,
        acceleration=acceleration,
        height=height)


def test_basic_pressure(test_args: Args) -> None:
    result = pressure_law.calculate_pressure(test_args.density_liquid,
        test_args.acceleration, test_args.height)
    assert_equal(result, 9827 * units.pascal)


def test_bad_density_liquid(test_args: Args) -> None:
    density_liquid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure(density_liquid, test_args.acceleration,
            test_args.height)
    with raises(TypeError):
        pressure_law.calculate_pressure(100, test_args.acceleration,
            test_args.height)


def test_bad_acceleration(test_args: Args) -> None:
    acceleration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure(test_args.density_liquid, acceleration,
            test_args.height)
    with raises(TypeError):
        pressure_law.calculate_pressure(test_args.density_liquid, 100,
            test_args.height)


def test_bad_height(test_args: Args) -> None:
    height = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure(test_args.density_liquid,
            test_args.acceleration, height)
    with raises(TypeError):
        pressure_law.calculate_pressure(test_args.density_liquid,
            test_args.acceleration, 100)
