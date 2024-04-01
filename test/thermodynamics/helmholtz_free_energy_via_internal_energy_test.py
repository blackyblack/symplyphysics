from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    helmholtz_free_energy_via_internal_energy as free_energy_law,)

# Description
## The internal energy of a system in thermodynamic equilibrium is 1 J, its temperature is
## 100 K and its entropy is 2e-2 J/K. Then its Helmholtz free energy is -1 J.

Args = namedtuple("Args", "u t s")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    u = Quantity(1 * units.joule)
    t = Quantity(100 * units.kelvin)
    s = Quantity(2e-2 * units.joule / units.kelvin)
    return Args(u=u, t=t, s=s)


def test_law(test_args: Args) -> None:
    result = free_energy_law.calculate_helmholtz_free_energy(test_args.u, test_args.t, test_args.s)
    assert_equal(result, -1 * units.joule)


def test_bad_energy(test_args: Args) -> None:
    ub = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        free_energy_law.calculate_helmholtz_free_energy(ub, test_args.t, test_args.s)
    with raises(TypeError):
        free_energy_law.calculate_helmholtz_free_energy(100, test_args.t, test_args.s)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        free_energy_law.calculate_helmholtz_free_energy(test_args.u, tb, test_args.s)
    with raises(TypeError):
        free_energy_law.calculate_helmholtz_free_energy(test_args.u, 100, test_args.s)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        free_energy_law.calculate_helmholtz_free_energy(test_args.u, test_args.t, sb)
    with raises(TypeError):
        free_energy_law.calculate_helmholtz_free_energy(test_args.u, test_args.t, 100)
