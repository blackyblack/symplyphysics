from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import kinetic_energy_via_momentum as kinetic_law

Args = namedtuple("Args", "p m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(1 * units.planck_momentum)
    m = Quantity(1 * units.planck_mass)
    return Args(p=p, m=m)


def test_law(test_args: Args) -> None:
    result = kinetic_law.calculate_kinetic_energy(test_args.p, test_args.m)
    assert_equal(result, 0.5 * units.planck_energy)


def test_bad_momentum(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        kinetic_law.calculate_kinetic_energy(pb, test_args.m)
    with raises(TypeError):
        kinetic_law.calculate_kinetic_energy(100, test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        kinetic_law.calculate_kinetic_energy(test_args.p, mb)
    with raises(TypeError):
        kinetic_law.calculate_kinetic_energy(test_args.p, 100)
