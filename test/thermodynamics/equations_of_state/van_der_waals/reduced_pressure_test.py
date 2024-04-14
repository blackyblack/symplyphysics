from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    prefixes,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import reduced_pressure

# Description
## The pressure of a van der Waals gas is 1 Pa, the value of critical pressure is 0.5 Pa. Then the reduced
## pressure of the gas is 2.

Args = namedtuple("Args", "p pc")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(1 * units.pascal)
    pc = Quantity(0.5 * units.pascal)
    return Args(p=p, pc=pc)


def test_law(test_args: Args) -> None:
    result = reduced_pressure.calculate_reduced_pressure(test_args.p, test_args.pc)
    assert_equal(result, 2)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reduced_pressure.calculate_reduced_pressure(pb, test_args.pc)
    with raises(TypeError):
        reduced_pressure.calculate_reduced_pressure(100, test_args.pc)
    with raises(errors.UnitsError):
        reduced_pressure.calculate_reduced_pressure(test_args.p, pb)
    with raises(TypeError):
        reduced_pressure.calculate_reduced_pressure(test_args.p, 100)
