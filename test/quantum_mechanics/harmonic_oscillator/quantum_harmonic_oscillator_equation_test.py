from collections import namedtuple
from pytest import fixture
from sympy import pi, exp, Rational, integrate, S
from symplyphysics import quantities
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.quantum_mechanics.harmonic_oscillator import quantum_harmonic_oscillator_equation

# Note
## The values are assumed to be in Planck units

Args = namedtuple("Args", "m w e psi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = 1
    w = 1
    e = Rational(3, 2)
    x = quantum_harmonic_oscillator_equation.position
    psi = (4 / pi)**Rational(1, 4) * x * exp(-1 * x**2 / 2)
    return Args(m=m, w=w, e=e, psi=psi)


def test_law(test_args: Args) -> None:
    # The wave function must be normalized
    assert expr_equals(
        integrate(abs(test_args.psi)**2, (quantum_harmonic_oscillator_equation.position, S.NegativeInfinity, S.Infinity)),
        1,
    )

    values = {
        quantities.hbar: 1,
        quantum_harmonic_oscillator_equation.particle_mass: test_args.m,
        quantum_harmonic_oscillator_equation.particle_energy: test_args.e,
        quantum_harmonic_oscillator_equation.angular_frequency: test_args.w,
        quantum_harmonic_oscillator_equation.wave_function(quantum_harmonic_oscillator_equation.position): test_args.psi,
    }
    lhs = quantum_harmonic_oscillator_equation.law.lhs.subs(values).doit()
    rhs = quantum_harmonic_oscillator_equation.law.rhs.subs(values).doit()
    assert expr_equals(lhs, rhs)
