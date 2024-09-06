from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.electricity import power_via_voltage_and_resistance as law

Args = namedtuple("Args", "v r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v = Quantity(4 * units.volt)
    r = Quantity(8 * units.ohm)
    return Args(v=v, r=r)


def test_law(test_args: Args) -> None:
    result = law.calculate_power(test_args.v, test_args.r)
    assert_equal(result, 2 * units.watt)


def test_bad_voltage(test_args: Args) -> None:
    vb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_power(vb, test_args.r)
    with raises(TypeError):
        law.calculate_power(100, test_args.r)


def test_bad_resistance(test_args: Args) -> None:
    rb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_power(test_args.v, rb)
    with raises(TypeError):
        law.calculate_power(test_args.v, 100)
