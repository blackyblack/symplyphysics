from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import mechanical_energy_is_kinetic_and_potential_energy as mechanical_energy_def

Args = namedtuple("Args", ["K", "P"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    K = Quantity(1 * units.joule)
    P = Quantity(5 * units.joule)
    return Args(K=K, P=P)


def test_basic_mechanical_energy(test_args: Args) -> None:
    result = mechanical_energy_def.calculate_mechanical_energy(test_args.K, test_args.P)
    assert_equal(result, 6 * units.joule)


def test_bad_energy(test_args: Args) -> None:
    Eb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        mechanical_energy_def.calculate_mechanical_energy(Eb, test_args.P)
    with raises(TypeError):
        mechanical_energy_def.calculate_mechanical_energy(100, test_args.P)
    with raises(errors.UnitsError):
        mechanical_energy_def.calculate_mechanical_energy(test_args.K, Eb)
    with raises(TypeError):
        mechanical_energy_def.calculate_mechanical_energy(test_args.K, 100)
