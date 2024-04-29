from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    quantum_isochoric_molar_heat_capacity_of_solids as einstein_law,
)

# Description
## The isochoric molar heat capacity for a solid in the Einstein model with atomic thermal oscillator
## frequency nu = 4 THz at temperature T = 300 K is C_V = 24.1 J/(K*mol).

Args = namedtuple("Args", "nu t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    nu = Quantity(4 * prefixes.tera * units.hertz)
    t = Quantity(300 * units.kelvin)
    return Args(nu=nu, t=t)


def test_law(test_args: Args) -> None:
    result = einstein_law.calculate_isochoric_molar_heat_capacity(test_args.nu, test_args.t)
    assert_equal(result, 24.1 * units.joule / (units.kelvin * units.mole))


def test_bad_frequency(test_args: Args) -> None:
    nub = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        einstein_law.calculate_isochoric_molar_heat_capacity(nub, test_args.t)
    with raises(TypeError):
        einstein_law.calculate_isochoric_molar_heat_capacity(100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        einstein_law.calculate_isochoric_molar_heat_capacity(test_args.nu, tb)
    with raises(TypeError):
        einstein_law.calculate_isochoric_molar_heat_capacity(test_args.nu, 100)
