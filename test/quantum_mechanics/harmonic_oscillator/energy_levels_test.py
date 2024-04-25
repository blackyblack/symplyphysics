from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
    prefixes,
)
from symplyphysics.laws.quantum_mechanics.harmonic_oscillator import energy_levels as energy_law

# Description
## The energy of the third level of a quantum harmonic oscillator for an angular frequency of
## 3.14 rad/fs is E_3 = 7.23 eV.

Args = namedtuple("Args", "n w")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n = 3
    w = Quantity(3.14 * units.radian / (prefixes.femto * units.second))
    return Args(n=n, w=w)


def test_law(test_args: Args) -> None:
    result = energy_law.calculate_energy_level(test_args.n, test_args.w)
    assert_equal(result, 7.23 * units.electronvolt)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy_level(nb, test_args.w)


def test_bad_angular_frequency(test_args: Args) -> None:
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy_level(test_args.n, wb)
    with raises(TypeError):
        energy_law.calculate_energy_level(test_args.n, 100)
