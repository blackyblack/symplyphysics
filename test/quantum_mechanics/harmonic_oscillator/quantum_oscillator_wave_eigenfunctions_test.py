from collections import namedtuple
from pytest import fixture, raises
from sympy import sqrt
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.core.symbols.quantities import evaluate_quantity
from symplyphysics.laws.quantum_mechanics.harmonic_oscillator import wave_eigenfunctions as wave_law

Args = namedtuple("Args", "n m w x")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n = 1

    # `sympy` evaluation halts unless these quantities have been evaluated to floating point form
    m = evaluate_quantity(units.planck_mass)
    w = evaluate_quantity(units.planck_angular_frequency)
    x = evaluate_quantity(9 * units.planck_length)
    return Args(n=n, m=m, w=w, x=x)


def test_law(test_args: Args) -> None:
    result = wave_law.calculate_wave_function_value(test_args.n, test_args.m, test_args.w, test_args.x)
    assert_equal(result, 6.13 / sqrt(units.meter))


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wave_law.calculate_wave_function_value(nb, test_args.m, test_args.w, test_args.x)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wave_law.calculate_wave_function_value(test_args.n, mb, test_args.w, test_args.x)
    with raises(TypeError):
        wave_law.calculate_wave_function_value(test_args.n, 100, test_args.w, test_args.x)


def test_bad_frequency(test_args: Args) -> None:
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wave_law.calculate_wave_function_value(test_args.n, test_args.m, wb, test_args.x)
    with raises(TypeError):
        wave_law.calculate_wave_function_value(test_args.n, test_args.m, 100, test_args.x)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wave_law.calculate_wave_function_value(test_args.n, test_args.m, test_args.w, xb)
    with raises(TypeError):
        wave_law.calculate_wave_function_value(test_args.n, test_args.m, test_args.w, 100)
