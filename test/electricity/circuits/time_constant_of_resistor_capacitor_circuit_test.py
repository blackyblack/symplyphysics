from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import time_constant_of_resistor_capacitor_circuit as law

Args = namedtuple("Args", "r c")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    r = Quantity(100 * units.ohm)
    c = Quantity(1e-6 * units.farad)
    return Args(r=r, c=c)


def test_law(test_args: Args) -> None:
    result = law.calculate_time_constant(test_args.r, test_args.c)
    assert_equal(result, 1e-4 * units.second)


def test_bad_resistance(test_args: Args) -> None:
    rb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_time_constant(rb, test_args.c)
    with raises(TypeError):
        law.calculate_time_constant(100, test_args.c)


def test_bad_conductance(test_args: Args) -> None:
    cb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_time_constant(test_args.r, cb)
    with raises(TypeError):
        law.calculate_time_constant(test_args.r, 100)
