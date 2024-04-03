from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.equations_of_state import virial_equation

# Description
## The second and third virial coefficients of a fluid are -800 cm**3/mol and -4 liter**2/mol**2 respectively.
## At molar density rho = 41 mol/m**3 the compressibility factor is Z = 0.960

Args = namedtuple("Args", "b c rho")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    b = Quantity(-800 * units.centimeter**3 / units.mole)
    c = Quantity(-4 * units.liter**2 / units.mole**2)
    rho = Quantity(41 * units.mole / units.meter**3)
    return Args(b=b, c=c, rho=rho)


def test_law(test_args: Args) -> None:
    result = virial_equation.calculate_compressibility_factor(test_args.b, test_args.c, test_args.rho)
    assert_equal(result, 0.96)


def test_bad_second_virial_coefficient(test_args: Args) -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        virial_equation.calculate_compressibility_factor(bb, test_args.c, test_args.rho)
    with raises(TypeError):
        virial_equation.calculate_compressibility_factor(100, test_args.c, test_args.rho)


def test_bad_third_virial_coefficient(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        virial_equation.calculate_compressibility_factor(test_args.b, cb, test_args.rho)
    with raises(TypeError):
        virial_equation.calculate_compressibility_factor(test_args.b, 100, test_args.rho)


def test_bad_molar_density(test_args: Args) -> None:
    rhob = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        virial_equation.calculate_compressibility_factor(test_args.b, test_args.c, rhob)
    with raises(TypeError):
        virial_equation.calculate_compressibility_factor(test_args.b, test_args.c, 100)
