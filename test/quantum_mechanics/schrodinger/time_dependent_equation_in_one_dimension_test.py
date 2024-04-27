from collections import namedtuple
from pytest import fixture
from sympy import exp, I, pi
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.quantum_mechanics.schrodinger import (
    time_dependent_equation_in_one_dimension as wave_function_eqn,
)

Args = namedtuple("Args", "m u psi x t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    position = wave_function_eqn.position
    time = wave_function_eqn.time
    m = 1
    u = 0
    psi = -3 * exp(I * (position - time / 2))
    x = 0.1
    t = 2
    return Args(m=m, u=u, psi=psi, x=x, t=t)


def test_law(test_args: Args) -> None:
    values = {
        wave_function_eqn.hbar: 1,
        wave_function_eqn.particle_mass: test_args.m,
        wave_function_eqn.potential_energy(wave_function_eqn.position, wave_function_eqn.time): test_args.u,
        wave_function_eqn.wave_function(wave_function_eqn.position, wave_function_eqn.time): test_args.psi,
    }
    lhs = wave_function_eqn.law.lhs.subs(values).doit()
    rhs = wave_function_eqn.law.rhs.subs(values).doit()
    assert expr_equals(lhs, rhs)
    values = {
        wave_function_eqn.position: test_args.x,
        wave_function_eqn.time: test_args.t
    }
    lhs = lhs.subs(values)
    rhs = rhs.subs(values)
    assert expr_equals(lhs, rhs)
