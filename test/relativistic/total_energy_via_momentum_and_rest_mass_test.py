from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.relativistic import total_energy_via_momentum_and_rest_mass as relation

# For a neutrino of rest mass m0 = 1.43e-36 kg having momentum p = 1.33e-27 kg*m/s
# the relativistic energy amounts to E = 2.61 eV.

Args = namedtuple("Args", "p m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(1.33e-27 * units.kilogram * units.meter / units.second)
    m = Quantity(1.43e-36 * units.kilogram)
    return Args(p=p, m=m)


def test_law(test_args: Args) -> None:
    result = relation.calculate_relativistic_energy(test_args.p, test_args.m)
    assert_equal(result, 2.61 * units.electronvolt, tolerance=2e-3)


def test_bad_momentum(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relation.calculate_relativistic_energy(pb, test_args.m)
    with raises(TypeError):
        relation.calculate_relativistic_energy(100, test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relation.calculate_relativistic_energy(test_args.p, mb)
    with raises(TypeError):
        relation.calculate_relativistic_energy(test_args.p, 100)
