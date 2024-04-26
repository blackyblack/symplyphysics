from collections import namedtuple
from pytest import fixture
from sympy import sin, pi, sqrt
from symplyphysics import units
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.quantum_mechanics.schrodinger import (
    time_independent_equation_in_one_dimension as wave_function_eqn,
)

Args = namedtuple("Args", "psi m u e")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    x = wave_function_eqn.position
    psi = sqrt(2) * sin(3 * pi * x)
    m = 1
    u = 1
    e = u + 9 * (pi * units.hbar)**2 / 2
    return Args(psi=psi, m=m, u=u, e=e)


def test_law(test_args: Args) -> None:
    values = {
        wave_function_eqn.wave_function(wave_function_eqn.position): test_args.psi,
        wave_function_eqn.potential_energy(wave_function_eqn.position): test_args.u,
        wave_function_eqn.particle_mass: test_args.m,
        wave_function_eqn.particle_energy: test_args.e
    }
    lhs = wave_function_eqn.law.lhs.subs(values)
    rhs = wave_function_eqn.law.rhs.subs(values)
    assert expr_equals(lhs, rhs)
