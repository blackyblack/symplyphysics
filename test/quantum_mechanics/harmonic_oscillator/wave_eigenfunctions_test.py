from collections import namedtuple
from pytest import fixture, raises
from sympy import sqrt
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
    prefixes,
)
from symplyphysics.core.symbols.quantities import evaluate_quantity
from symplyphysics.laws.quantum_mechanics.harmonic_oscillator import wave_eigenfunctions as wave_law

Args = namedtuple("Args", "n m w x")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n = 1
    m = evaluate_quantity(units.planck_mass)
    w = evaluate_quantity(units.planck_angular_frequency)
    x = evaluate_quantity(9 * units.planck_length)
    return Args(n=n, m=m, w=w, x=x)


def test_law(test_args: Args) -> None:
    result = wave_law.calculate_wave_function(test_args.n, test_args.m, test_args.w, test_args.x)
    assert_equal(result, 6.13 / sqrt(units.meter))
