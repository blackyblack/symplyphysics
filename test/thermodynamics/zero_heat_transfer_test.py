from collections import namedtuple
from pytest import fixture, raises
from sympy import Rational
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import zero_heat_transfer

Args = namedtuple("Args", ["n", "t0", "V0", "V1", "y"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n = Quantity(1 * units.mole)
    t0 = Quantity(100 * units.kelvin)
    V0 = Quantity(1000 * units.liter)
    V1 = Quantity(2000 * units.liter)
    # Choose specific heats ratio
    y = Rational(1.665)
    return Args(n=n, t0=t0, V0=V0, V1=V1, y=y)


def test_basic_pressure(test_args: Args) -> None:
    result = zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, test_args.V0,
        test_args.V1, test_args.y)
    assert_equal(result, 262.19 * units.pascal)


def test_bad_mole_count(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(nb, test_args.t0, test_args.V0, test_args.V1,
            test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(100, test_args.t0, test_args.V0, test_args.V1,
            test_args.y)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(test_args.n, tb, test_args.V0, test_args.V1,
            test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(test_args.n, 100, test_args.V0, test_args.V1,
            test_args.y)


def test_bad_volume(test_args: Args) -> None:
    Vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, Vb, test_args.V1,
            test_args.y)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, test_args.V0, Vb,
            test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, 100, test_args.V1,
            test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, test_args.V0, 100,
            test_args.y)


def test_bad_specific_heats_ratio(test_args: Args) -> None:
    yb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, test_args.V0, test_args.V1,
            yb)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, test_args.V0, test_args.V1,
            'bad')
