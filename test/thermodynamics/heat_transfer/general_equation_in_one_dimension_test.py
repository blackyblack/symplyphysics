from collections import namedtuple
from pytest import fixture
from sympy import sin, exp
from symplyphysics import units, Quantity
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.heat_transfer import (
    general_equation_in_one_dimension as heat_equation
)

Args = namedtuple("Args", "rho cv k q T t x")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    rho = 1
    cv = 1
    k = 1
    q = 0
    T = sin(heat_equation.position) * exp(-1 * heat_equation.time)
    t = 1
    x = 1
    return Args(rho=rho, cv=cv, k=k, q=q, T=T, t=t, x=x)


def test_law(test_args: Args) -> None:
    values = {
        heat_equation.medium_density: test_args.rho,
        heat_equation.medium_specific_heat_capacity: test_args.cv,
        heat_equation.thermal_conductivity(heat_equation.position): test_args.k,
        heat_equation.heat_source_density(heat_equation.position): test_args.q,
        heat_equation.temperature(heat_equation.position, heat_equation.time): test_args.T,
    }
    lhs = heat_equation.law.lhs.subs(values)
    rhs = heat_equation.law.rhs.subs(values)
    assert expr_equals(lhs, rhs)
    values = {
        heat_equation.position: test_args.x,
        heat_equation.time: test_args.t,
    }
    lhs = lhs.subs(values)
    rhs = rhs.subs(values)
    assert expr_equals(lhs, rhs)
