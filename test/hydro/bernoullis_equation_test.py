from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
import symplyphysics.laws.hydro.bernoullis_equation as bernoullis_equation

# Description
# rho = 1 kg/m^3
# Initial parameters: P = 1 Pa, v = 2 m/s, h = 0 m
# Final parameters: P = 1.5 Pa, v = 1 m/s, h = 0.102 m
# The Bernoulli's equation should verify these figures.


@fixture(name="test_args")
def test_args_fixture():
    P1 = Quantity(1 * units.pascal)
    rho = Quantity(1 * units.kilogram / units.meter**3)
    v1 = Quantity(2 * units.meter / units.second)
    h1 = Quantity(0 * units.meter)
    v2 = Quantity(1 * units.meter / units.second)
    h2 = Quantity(0.102 * units.meter)
    Args = namedtuple("Args", ["P1", "rho", "v1", "h1", "v2", "h2"])
    return Args(P1=P1, rho=rho, v1=v1, h1=h1, v2=v2, h2=h2)


def test_bernoullis_equation(test_args):
    result = bernoullis_equation.calculate_pressure(
        test_args.P1, test_args.rho, test_args.v1, test_args.h1, test_args.v2, test_args.h2
    )
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_pressure = convert_to(result, units.pascal).evalf(3)
    assert result_pressure == approx(1.5, 0.003)


def test_bad_pressure_before(test_args):
    P1_bad = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bernoullis_equation.calculate_pressure(
            P1_bad, test_args.rho, test_args.v1, test_args.h1, test_args.v2, test_args.h2
        )
    with raises(TypeError):
        bernoullis_equation.calculate_pressure(
            100, test_args.rho, test_args.v1, test_args.h1, test_args.v2, test_args.h2
        )


def test_bad_density(test_args):
    rho_bad = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bernoullis_equation.calculate_pressure(
            test_args.P1, rho_bad, test_args.v1, test_args.h1, test_args.v2, test_args.h2
        )
    with raises(TypeError):
        bernoullis_equation.calculate_pressure(
            test_args.P1, 100, test_args.v1, test_args.h1, test_args.v2, test_args.h2
        )


def test_bad_speed_before(test_args):
    v1_bad = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bernoullis_equation.calculate_pressure(
            test_args.P1, test_args.rho, v1_bad, test_args.h1, test_args.v2, test_args.h2
        )
    with raises(TypeError):
        bernoullis_equation.calculate_pressure(
            test_args.P1, test_args.rho, 100, test_args.h1, test_args.v2, test_args.h2
        )


def test_bad_elevation_before(test_args):
    h1_bad = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bernoullis_equation.calculate_pressure(
            test_args.P1, test_args.rho, test_args.v1, h1_bad, test_args.v2, test_args.h2
        )
    with raises(TypeError):
        bernoullis_equation.calculate_pressure(
            test_args.P1, test_args.rho, test_args.v1, 100, test_args.v2, test_args.h2
        )


def test_bad_speed_after(test_args):
    v2_bad = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bernoullis_equation.calculate_pressure(
            test_args.P1, test_args.rho, test_args.v1, test_args.h1, v2_bad, test_args.h2
        )
    with raises(TypeError):
        bernoullis_equation.calculate_pressure(
            test_args.P1, test_args.rho, test_args.v1, test_args.h1, 100, test_args.h2
        )


def test_bad_elevation_after(test_args):
    h2_bad = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bernoullis_equation.calculate_pressure(
            test_args.P1, test_args.rho, test_args.v1, test_args.h1, test_args.v2, h2_bad
        )
    with raises(TypeError):
        bernoullis_equation.calculate_pressure(
            test_args.P1, test_args.rho, test_args.v1, test_args.h1, test_args.v2, 100
        )
