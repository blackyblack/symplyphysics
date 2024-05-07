from collections import namedtuple
from pytest import fixture, raises
from sympy import sqrt, I
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.quantum_mechanics import probability_density_of_quantum_state as probability_law

Args = namedtuple("Args", "psi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    psi = Quantity((1 + 3*I) / sqrt(units.meter))
    return Args(psi=psi)


def test_law(test_args: Args) -> None:
    result = probability_law.calculate_probability_density(test_args.psi)
    assert_equal(result, 10 / units.meter)


def test_bad_wave_function() -> None:
    psib = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        probability_law.calculate_probability_density(psib)
    with raises(TypeError):
        probability_law.calculate_probability_density(100)
