from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity.radial_motion import total_energy_is_negative_average_kinetic_energy as law

Args = namedtuple("Args", "k")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k = Quantity(3e8 * units.joule)
    return Args(k=k)


def test_law(test_args: Args) -> None:
    result = law.calculate_total_energy(test_args.k)
    assert_equal(result, -3e8 * units.joule)


def test_bad_energy() -> None:
    kb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_total_energy(kb)
    with raises(TypeError):
        law.calculate_total_energy(100)
