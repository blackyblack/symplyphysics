from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity.radial_motion import potential_energy_of_planetary_motion as law

Args = namedtuple("Args", "u l m r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    u = Quantity(-5.3e33 * units.joule)
    l = Quantity(7.8e33 * units.kilogram * units.meter**2 / units.second)
    m = Quantity(5.97e33 * units.kilogram)
    r = Quantity(1.0 * units.astronomical_unit)
    return Args(u=u, l=l, m=m, r=r)


def test_law(test_args: Args) -> None:
    result = law.calculate_planetary_potential_energy(test_args.u, test_args.l, test_args.m, test_args.r)
    assert_equal(result, -5.3e33 * units.joule)


def test_bad_energy(test_args: Args) -> None:
    ub = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_planetary_potential_energy(ub, test_args.l, test_args.m, test_args.r)
    with raises(TypeError):
        law.calculate_planetary_potential_energy(100, test_args.l, test_args.m, test_args.r)


def test_bad_angular_momentum(test_args: Args) -> None:
    lb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_planetary_potential_energy(test_args.u, lb, test_args.m, test_args.r)
    with raises(TypeError):
        law.calculate_planetary_potential_energy(test_args.u, 100, test_args.m, test_args.r)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_planetary_potential_energy(test_args.u, test_args.l, mb, test_args.r)
    with raises(TypeError):
        law.calculate_planetary_potential_energy(test_args.u, test_args.l, 100, test_args.r)


def test_bad_distance(test_args: Args) -> None:
    rb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_planetary_potential_energy(test_args.u, test_args.l, test_args.m, rb)
    with raises(TypeError):
        law.calculate_planetary_potential_energy(test_args.u, test_args.l, test_args.m, 100)
