from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import laplase_pressure as laplase_law

# Description
## Test example from http://ru.solverbook.com/spravochnik/molekulyarnaya-fizika-i-termodinamika/davlenie-pod-iskrivlennoj-poverxnostyu-zhidkosti/
## Water surface tension coefficient is equal 7.4 * 10^(-2) N/m. The radius of the bubble in the water is equal 0,005 mm. Then the Laplace pressure will equal 2.96 * 10^4 Pa.

Args = namedtuple("Args", ["sigma", "r"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    sigma = Quantity(7.4e-2 * units.newton / units.meter)
    r = Quantity(5e-3 * units.millimeters)
    return Args(sigma=sigma, r=r)


def test_basic_laplase_pressure(test_args: Args) -> None:
    result = laplase_law.calculate_laplase_pressure(test_args.sigma, test_args.r)
    assert_equal(result, 2.96e4 * units.pascal)


def test_bad_surface_tension(test_args: Args) -> None:
    sigma_b = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        laplase_law.calculate_laplase_pressure(sigma_b, test_args.r)
    with raises(TypeError):
        laplase_law.calculate_laplase_pressure(100, test_args.r)


def test_bad_radius(test_args: Args) -> None:
    rb = Quantity(10 * units.coulomb)
    with raises(errors.UnitsError):
        laplase_law.calculate_laplase_pressure(test_args.sigma, rb)
    with raises(TypeError):
        laplase_law.calculate_laplase_pressure(test_args.sigma, 100)
