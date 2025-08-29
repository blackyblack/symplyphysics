from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.electricity import current_is_voltage_over_impedance as law

Args = namedtuple("Args", ["v", "z"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v = Quantity((3 - 1j) * units.volt)
    z = Quantity((4 + 3j) * units.ohm)
    return Args(v=v, z=z)


def test_basic_current(test_args: Args) -> None:
    result = law.calculate_current(test_args.v, test_args.z)
    assert_equal(result, (0.36 - 0.52j) * units.ampere)


def test_bad_voltage(test_args: Args) -> None:
    vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        law.calculate_current(vb, test_args.z)
    with raises(TypeError):
        law.calculate_current(100, test_args.z)


def test_bad_impedance(test_args: Args) -> None:
    zb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        law.calculate_current(test_args.v, zb)
    with raises(TypeError):
        law.calculate_current(test_args.v, 100)
