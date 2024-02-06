from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, Quantity, units, errors
from symplyphysics.laws.hydro import reynolds_number

# Description
# There is pipe with diameter 0.1 m, density 1000 kg/m3, velocity 1 m/s
# and dynamic viscosity 0.000894 Pa*s. The reynolds number should be 111856.823


@fixture(name="test_args")
def test_args_fixture():
    rho = Quantity(1000 * units.kilogram / units.meter**3)
    d = Quantity(0.1 * units.meter)
    v = Quantity(1 * units.meter / units.second)
    mu = Quantity(0.000894 * units.pascal * units.second)
    Args = namedtuple("Args", ["d", "rho", "v", "mu"])
    return Args(d=d, rho=rho, v=v, mu=mu)


def test_reynolds_number(test_args):
    result = reynolds_number.calculate_reynolds_number(test_args.d, test_args.rho, test_args.v,
        test_args.mu)
    assert_equal(result, 111856.823)


def test_bad_velocity(test_args):
    bv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reynolds_number.calculate_reynolds_number(test_args.d, test_args.rho, bv, test_args.mu)
    with raises(TypeError):
        reynolds_number.calculate_reynolds_number(test_args.d, test_args.rho, 0, test_args.mu)


def test_bad_diameter(test_args):
    bd = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reynolds_number.calculate_reynolds_number(bd, test_args.rho, test_args.v, test_args.mu)
    with raises(TypeError):
        reynolds_number.calculate_reynolds_number(0, test_args.rho, test_args.v, test_args.mu)


def test_bad_density(test_args):
    bd = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reynolds_number.calculate_reynolds_number(test_args.d, bd, test_args.v, test_args.mu)
    with raises(TypeError):
        reynolds_number.calculate_reynolds_number(test_args.d, 10, test_args.v, test_args.mu)


def test_bad_dynamic_viscosity(test_args):
    bd = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reynolds_number.calculate_reynolds_number(test_args.d, test_args.rho, test_args.v, bd)
    with raises(TypeError):
        reynolds_number.calculate_reynolds_number(test_args.d, test_args.rho, test_args.v, 0.05)
