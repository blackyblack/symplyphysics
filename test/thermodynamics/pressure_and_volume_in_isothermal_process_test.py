from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import pressure_and_volume_in_isothermal_process as boyles_law

Args = namedtuple("Args", ["P0", "P1", "V0"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    P0 = Quantity(1 * units.pascal)
    P1 = Quantity(2 * units.pascal)
    V0 = Quantity(1 * units.liter)
    return Args(P0=P0, P1=P1, V0=V0)


def test_basic_volume(test_args: Args) -> None:
    result = boyles_law.calculate_volume(test_args.P0, test_args.V0, test_args.P1)
    assert_equal(result, 0.5 * units.liter)


def test_bad_pressure(test_args: Args) -> None:
    Pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        boyles_law.calculate_volume(Pb, test_args.V0, test_args.P1)
    with raises(errors.UnitsError):
        boyles_law.calculate_volume(test_args.P0, test_args.V0, Pb)
    with raises(TypeError):
        boyles_law.calculate_volume(100, test_args.V0, test_args.P1)
    with raises(TypeError):
        boyles_law.calculate_volume(test_args.P0, test_args.V0, 100)


def test_bad_volume(test_args: Args) -> None:
    Vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        boyles_law.calculate_volume(test_args.P0, Vb, test_args.P1)
    with raises(TypeError):
        boyles_law.calculate_volume(test_args.P0, 100, test_args.P1)
