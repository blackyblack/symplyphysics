from collections import namedtuple
from pytest import fixture
from sympy import exp, I, symbols, sqrt
from symplyphysics import quantities
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.quantum_mechanics.schrodinger import (
    time_dependent_equation_in_one_dimension as wave_function_eqn,)

Args = namedtuple("Args", "u psi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    x = wave_function_eqn.position
    t = wave_function_eqn.time
    m = wave_function_eqn.particle_mass
    e, u0, w = symbols("E U_0 Omega", positive=True)
    u = u0 * exp(I * w * t)
    phase = sqrt(2 * m * e) * x - e * t + I * (u0 / w) * exp(I * w * t)
    psi = exp((I / quantities.hbar) * phase)
    return Args(u=u, psi=psi)


def test_law(test_args: Args) -> None:
    values = {
        wave_function_eqn.potential_energy(wave_function_eqn.position, wave_function_eqn.time):
        test_args.u,
        wave_function_eqn.wave_function(wave_function_eqn.position, wave_function_eqn.time):
        test_args.psi,
    }
    lhs = wave_function_eqn.law.lhs.subs(values).doit()
    rhs = wave_function_eqn.law.rhs.subs(values).doit()
    assert expr_equals(lhs, rhs)
