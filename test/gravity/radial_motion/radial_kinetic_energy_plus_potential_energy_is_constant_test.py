from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity.radial_motion import (
    radial_kinetic_energy_plus_potential_energy_is_constant as law,)

Args = namedtuple("Args", "m v u")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(5.972e24 * units.kilogram)
    v = Quantity(30 * units.kilometer / units.second)
    u = Quantity(-5.3e33 * units.joule)
    return Args(m=m, v=v, u=u)


def test_law(test_args: Args) -> None:
    result = law.calculate_total_energy(test_args.m, test_args.v, test_args.u)
    assert_equal(result, -2.6e33 * units.joule, tolerance=5e-3)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_total_energy(mb, test_args.v, test_args.u)
    with raises(TypeError):
        law.calculate_total_energy(100, test_args.v, test_args.u)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_total_energy(test_args.m, vb, test_args.u)
    with raises(TypeError):
        law.calculate_total_energy(test_args.m, 100, test_args.u)


def test_bad_energy(test_args: Args) -> None:
    ub = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_total_energy(test_args.m, test_args.v, ub)
    with raises(TypeError):
        law.calculate_total_energy(test_args.m, test_args.v, 100)
