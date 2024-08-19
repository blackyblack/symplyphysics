from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import photon_energy_is_proportional_to_angular_frequency as planck_law

Args = namedtuple("Args", ["w"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w = Quantity(2 * pi * 3e16 * units.radian / units.second)
    return Args(w=w)


def test_basic_energy(test_args: Args) -> None:
    result = planck_law.calculate_energy(test_args.w)
    assert_equal(result, 1.988e-17 * units.joule)


def test_bad_frequency() -> None:
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        planck_law.calculate_energy(fb)
    with raises(TypeError):
        planck_law.calculate_energy(100)
