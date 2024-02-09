from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import volume_is_constant as isochoric_law

Args = namedtuple("Args", ["t0", "t1", "P0"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t0 = Quantity(1 * units.kelvin)
    t1 = Quantity(2 * units.kelvin)
    P0 = Quantity(1 * units.pascal)
    return Args(t0=t0, t1=t1, P0=P0)


def test_basic_pressure(test_args: Args) -> None:
    result = isochoric_law.calculate_pressure(test_args.t0, test_args.P0, test_args.t1)
    assert_equal(result, 2 * units.pascal)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        isochoric_law.calculate_pressure(tb, test_args.P0, test_args.t1)
    with raises(errors.UnitsError):
        isochoric_law.calculate_pressure(test_args.t0, test_args.P0, tb)
    with raises(TypeError):
        isochoric_law.calculate_pressure(100, test_args.P0, test_args.t1)
    with raises(TypeError):
        isochoric_law.calculate_pressure(test_args.t0, test_args.P0, 100)


def test_bad_pressure(test_args: Args) -> None:
    Pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        isochoric_law.calculate_pressure(test_args.t0, Pb, test_args.t1)
    with raises(TypeError):
        isochoric_law.calculate_pressure(test_args.t0, 100, test_args.t1)
