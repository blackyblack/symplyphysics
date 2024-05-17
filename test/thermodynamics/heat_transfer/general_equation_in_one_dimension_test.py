from collections import namedtuple
from pytest import fixture
from sympy import sin, exp, Rational
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.heat_transfer import (general_equation_in_one_dimension as
    heat_equation)

# Description
## The values of the parameters are that of carbon monoxide (CO).

Args = namedtuple("Args", "rho cv k q T")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    rho = Rational(15, 10)  # kg/m**3
    cv = 1500  # J/(K * kg)
    k = Rational(24, 1000)  # W/(K*m)
    q = 0
    T = sin(100 * heat_equation.position) * exp(Rational(-8, 75) * heat_equation.time)
    return Args(rho=rho, cv=cv, k=k, q=q, T=T)


def test_law(test_args: Args) -> None:
    values = {
        heat_equation.medium_density: test_args.rho,
        heat_equation.medium_specific_heat_capacity: test_args.cv,
        heat_equation.thermal_conductivity(heat_equation.position): test_args.k,
        heat_equation.heat_source_density(heat_equation.position, heat_equation.time): test_args.q,
        heat_equation.temperature(heat_equation.position, heat_equation.time): test_args.T,
    }
    lhs = heat_equation.law.lhs.subs(values).doit()
    rhs = heat_equation.law.rhs.subs(values).doit()
    assert expr_equals(lhs, rhs)
