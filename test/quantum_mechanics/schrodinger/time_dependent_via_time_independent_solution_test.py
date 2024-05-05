from collections import namedtuple
from pytest import fixture, raises
from sympy import sqrt
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.quantum_mechanics.schrodinger import (
    time_dependent_via_time_independent_solution as law
)

Args = namedtuple("Args", "psi e t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    psi = Quantity((1+3j) / sqrt(units.meter))
    e = Quantity(1 * units.electronvolt)
    t = Quantity(1e-3 * units.second)
    return Args(psi=psi, e=e, t=t)


def test_law(test_args: Args) -> None:
    result = law.calculate_time_dependent_wave_function_value(test_args.psi, test_args.e, test_args.t)
    assert_equal(result, (-0.844-3.05j) / sqrt(units.meter))


def test_bad_wave_function(test_args: Args) -> None:
    psib = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_time_dependent_wave_function_value(psib, test_args.e, test_args.t)
    with raises(TypeError):
        law.calculate_time_dependent_wave_function_value(100, test_args.e, test_args.t)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_time_dependent_wave_function_value(test_args.psi, eb, test_args.t)
    with raises(TypeError):
        law.calculate_time_dependent_wave_function_value(test_args.psi, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_time_dependent_wave_function_value(test_args.psi, test_args.e, tb)
    with raises(TypeError):
        law.calculate_time_dependent_wave_function_value(test_args.psi, test_args.e, 100)
