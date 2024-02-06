from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, Quantity, units, errors
from symplyphysics.laws.hydro import hagen_poiseuille_equation


@fixture(name="test_args")
def test_args_fixture():
    mu = Quantity(0.000894 * units.pascal * units.second)
    l = Quantity(0.1 * units.meter)
    q = Quantity(1 * (units.meter**3 / units.second))
    r = Quantity(0.1 * units.meter)
    Args = namedtuple("Args", ["mu", "l", "q", "r"])
    return Args(mu=mu, l=l, q=q, r=r)


def test_basic(test_args):
    result = hagen_poiseuille_equation.calculate_delta_pressure(test_args.mu, test_args.l,
        test_args.q, test_args.r)
    assert_equal(result, 2.276 * units.pascal)


def test_bad_dynamic_viscosity(test_args):
    bad_mu = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hagen_poiseuille_equation.calculate_delta_pressure(bad_mu, test_args.l, test_args.q,
            test_args.r)
    with raises(TypeError):
        hagen_poiseuille_equation.calculate_delta_pressure(0, test_args.l, test_args.q, test_args.r)


def test_bad_length(test_args):
    bad_l = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hagen_poiseuille_equation.calculate_delta_pressure(test_args.mu, bad_l, test_args.q,
            test_args.r)
    with raises(TypeError):
        hagen_poiseuille_equation.calculate_delta_pressure(test_args.mu, 0, test_args.q,
            test_args.r)


def test_bad_flow_rate(test_args):
    bad_q = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hagen_poiseuille_equation.calculate_delta_pressure(test_args.mu, test_args.l, bad_q,
            test_args.r)
    with raises(TypeError):
        hagen_poiseuille_equation.calculate_delta_pressure(test_args.mu, test_args.l, 0,
            test_args.r)


def test_bad_radius(test_args):
    bad_r = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hagen_poiseuille_equation.calculate_delta_pressure(test_args.mu, test_args.l, test_args.q,
            bad_r)
    with raises(TypeError):
        hagen_poiseuille_equation.calculate_delta_pressure(test_args.mu, test_args.l, test_args.q,
            0)
