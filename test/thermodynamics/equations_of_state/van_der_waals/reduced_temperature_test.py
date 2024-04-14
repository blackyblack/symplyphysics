from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import reduced_temperature

# Description
## The temperature of a van der Waals gas is 400 K, the value of critical temperature is 200 K. Then the reduced
## temperature of the gas is 2.

Args = namedtuple("Args", "t tc")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(400 * units.kelvin)
    tc = Quantity(200 * units.kelvin)
    return Args(t=t, tc=tc)


def test_law(test_args: Args) -> None:
    result = reduced_temperature.calculate_reduced_temperature(test_args.t, test_args.tc)
    assert_equal(result, 2)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reduced_temperature.calculate_reduced_temperature(tb, test_args.tc)
    with raises(TypeError):
        reduced_temperature.calculate_reduced_temperature(100, test_args.tc)
    with raises(errors.UnitsError):
        reduced_temperature.calculate_reduced_temperature(test_args.t, tb)
    with raises(TypeError):
        reduced_temperature.calculate_reduced_temperature(test_args.t, 100)
