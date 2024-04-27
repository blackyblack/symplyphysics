from collections import namedtuple
from pytest import fixture
from sympy import sin, exp
from symplyphysics import units, Quantity
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.heat_transfer import (
    equation_in_homogeneous_medium_in_one_dimension as heat_equation
)

Args = namedtuple("Args", "temperature chi x t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    chi = Quantity(95 * units.millimeter**2 / units.second)
    k = Quantity(0.5 / units.meter)
    temperature = sin(k * heat_equation.position) * exp(-chi * k**2 * heat_equation.time)
    x = Quantity(1 * units.meter)
    t = Quantity(0.1 * units.second)
    return Args(temperature=temperature, chi=chi, x=x, t=t)


def test_law(test_args: Args) -> None:
    values = {
        heat_equation.temperature(heat_equation.position, heat_equation.time): test_args.temperature,
        heat_equation.thermal_diffusivity: test_args.chi
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
