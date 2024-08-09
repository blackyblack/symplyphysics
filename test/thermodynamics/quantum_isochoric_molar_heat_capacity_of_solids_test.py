from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    quantum_isochoric_molar_heat_capacity_of_solids as einstein_law,)

# Description
## The isochoric molar heat capacity for a solid in the Einstein model with atomic thermal oscillator
## frequency nu = 4 THz at temperature T = 300 K is C_V = 24.1 J/(K*mol). The reduced photon energy
## is x = 0.64.

Args = namedtuple("Args", "x")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    nu = Quantity(4 * prefixes.tera * units.hertz)
    w = Quantity(2 * pi * nu)
    t = Quantity(300 * units.kelvin)
    x = Quantity(units.hbar * w / (units.boltzmann_constant * t))
    return Args(x=x)


def test_law(test_args: Args) -> None:
    result = einstein_law.calculate_isochoric_molar_heat_capacity(test_args.x)
    assert_equal(result, 24.1 * units.joule / (units.kelvin * units.mole))


def test_bad_number() -> None:
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        einstein_law.calculate_isochoric_molar_heat_capacity(xb)
