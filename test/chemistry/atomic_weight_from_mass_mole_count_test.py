from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)
from symplyphysics.laws.chemistry import atomic_weight_from_mass_mole_count

Args = namedtuple("Args", ["m", "N"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # molar mass of water is 18.0153 gram / mole
    mass = Quantity(18.0153 * units.gram)
    mole_count = Quantity(1 * units.mole)
    return Args(m=mass, N=mole_count)


def test_basic_atomic_weight(test_args: Args) -> None:
    result = atomic_weight_from_mass_mole_count.calculate_atomic_weight(test_args.m, test_args.N)
    assert_equal(result, 18.0153 * units.gram / units.mole)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(mb, test_args.N)
    with raises(TypeError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(100, test_args.N)


def test_bad_mole_count(test_args: Args) -> None:
    Nb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(test_args.m, Nb)
    with raises(TypeError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(test_args.m, 100)
