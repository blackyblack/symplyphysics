from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import electromotive_force_induced_in_rotating_coil as law

Args = namedtuple("Args", "n dphi dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n = 100
    dphi = Quantity(-5e-3 * units.weber)
    dt = Quantity(1e-4 * units.second)
    return Args(n=n, dphi=dphi, dt=dt)


def test_law(test_args: Args) -> None:
    result = law.calculate_electromotive_force(test_args.n, test_args.dphi, test_args.dt)
    assert_equal(result, 5e3 * units.volt)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_electromotive_force(nb, test_args.dphi, test_args.dt)


def test_bad_magnetic_flux(test_args: Args) -> None:
    phi_b = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_electromotive_force(test_args.n, phi_b, test_args.dt)
    with raises(TypeError):
        law.calculate_electromotive_force(test_args.n, 100, test_args.dt)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_electromotive_force(test_args.n, test_args.dphi, tb)
    with raises(TypeError):
        law.calculate_electromotive_force(test_args.n, test_args.dphi, 100)
