from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.thermodynamics import laplas_pressure as laplas_law

# Description
## Test example from http://ru.solverbook.com/spravochnik/molekulyarnaya-fizika-i-termodinamika/davlenie-pod-iskrivlennoj-poverxnostyu-zhidkosti/
## Water surface tension coefficient is equal 7.4 * 10^(-2) N/m. The radius of the bubble in the water is equal 0,005 mm. Then the Laplace pressure will equal 2.96 * 10^4 Pa.


@fixture(name="test_args")
def test_args_fixture():
    sigma = Quantity(7.4e-2 * units.newton / units.meter)
    r = Quantity(5e-3 * units.millimeters)
    Args = namedtuple("Args", ["sigma", "r"])
    return Args(sigma=sigma, r=r)


def test_basic_laplas_pressure(test_args):
    result = laplas_law.calculate_laplas_pressure(test_args.sigma, test_args.r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_laplas_pressure = convert_to(result, units.pascal).evalf(5)
    assert_approx(result_laplas_pressure, 2.96e4)


def test_bad_surface_tension(test_args):
    sigma_b = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        laplas_law.calculate_laplas_pressure(sigma_b, test_args.r)
    with raises(TypeError):
        laplas_law.calculate_laplas_pressure(100, test_args.r)


def test_bad_radius(test_args):
    rb = Quantity(10 * units.coulomb)
    with raises(errors.UnitsError):
        laplas_law.calculate_laplas_pressure(test_args.sigma, rb)
    with raises(TypeError):
        laplas_law.calculate_laplas_pressure(test_args.sigma, 100)
