from collections import namedtuple
from pytest import fixture
from sympy import sin, exp
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.heat_transfer import (
    equation_in_homogeneous_medium_in_one_dimension as heat_equation
)

Args = namedtuple("Args", "temperature chi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    chi = 1e-5
    temperature = sin(heat_equation.position) * exp(-chi * heat_equation.time)
    return Args(temperature=temperature, chi=chi)


def test_law(test_args: Args) -> None:
    values = {
        heat_equation.temperature(heat_equation.position, heat_equation.time): test_args.temperature,
        heat_equation.thermal_diffusivity: test_args.chi
    }
    lhs = heat_equation.law.lhs.subs(values)
    rhs = heat_equation.law.rhs.subs(values)
    assert expr_equals(lhs, rhs)
