from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import average_square_speed_in_maxwell_boltzmann_statistics as velocity_law

# Description
## Example from https://easyfizika.ru/zadachi/molekulyarnaya-fizika/opredelit-srednyuyu-kvadratichnuyu-skorost-molekul-vodoroda/
## The ideal gas has a temperature of 273 K. The mass of molecule is equal to 3,32 * 10^(-27) kg. Then the average square of velocity should be 3.403 * 10^6 (m/s)^2

Args = namedtuple("Args", ["m", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(3.32e-27 * units.kilogram)
    t = Quantity(273 * units.kelvin)
    return Args(m=m, t=t)


def test_basic_average_square_velocity(test_args: Args) -> None:
    result = velocity_law.calculate_average_square_velocity(test_args.t, test_args.m)
    assert_equal(result, 3.403e6 * (units.meter / units.second)**2)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_average_square_velocity(test_args.t, mb)
    with raises(TypeError):
        velocity_law.calculate_average_square_velocity(test_args.t, 100)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(100 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_average_square_velocity(tb, test_args.m)
    with raises(TypeError):
        velocity_law.calculate_average_square_velocity(100, test_args.m)
