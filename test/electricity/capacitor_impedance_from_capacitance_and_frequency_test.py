from collections import namedtuple
from sympy import I
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.electricity import capacitor_impedance_from_capacitance_and_frequency as capacitor_impedance_law

# Description
## Assert we have a capacitor with 5F capacitance. In 0.05Hz (0.314159 rad/s) circuit impedance of this element should be 636.61mOhm.
## (https://www.translatorscafe.com/unit-converter/ru-RU/calculator/capacitor-impedance/)

Args = namedtuple("Args", ["c", "w"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    c = Quantity(5 * units.farad)
    w = Quantity(0.314159 * units.radian / units.second)
    return Args(c=c, w=w)


def test_basic_impedance(test_args: Args) -> None:
    result = capacitor_impedance_law.calculate_impedance(test_args.c, test_args.w)
    assert_equal(result, -0.63662 * I * units.ohm)


def test_bad_capacitance(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitor_impedance_law.calculate_impedance(cb, test_args.w)
    with raises(TypeError):
        capacitor_impedance_law.calculate_impedance(100, test_args.w)


def test_bad_frequency(test_args: Args) -> None:
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitor_impedance_law.calculate_impedance(test_args.c, wb)
    with raises(TypeError):
        capacitor_impedance_law.calculate_impedance(test_args.c, 100)
