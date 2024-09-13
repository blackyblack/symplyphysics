from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import inductance_is_magnetic_flux_over_current as law

Args = namedtuple("Args", "phi i")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    phi = Quantity(4 * units.weber)
    i = Quantity(0.2 * units.ampere)
    return Args(phi=phi, i=i)


def test_law(test_args: Args) -> None:
    result = law.calculate_inductance(test_args.phi, test_args.i)
    assert_equal(result, 20 * units.henry)


def test_bad_magnetic_flux(test_args: Args) -> None:
    phi_b = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_inductance(phi_b, test_args.i)
    with raises(TypeError):
        law.calculate_inductance(100, test_args.i)


def test_bad_current(test_args: Args) -> None:
    ib = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_inductance(test_args.phi, ib)
    with raises(TypeError):
        law.calculate_inductance(test_args.phi, 100)
