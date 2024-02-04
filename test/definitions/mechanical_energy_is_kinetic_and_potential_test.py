from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import mechanical_energy_is_kinetic_and_potential as mechanical_energy_def


@fixture(name="test_args")
def test_args_fixture():
    K = Quantity(1 * units.joule)
    P = Quantity(5 * units.joule)
    Args = namedtuple("Args", ["K", "P"])
    return Args(K=K, P=P)


def test_basic_mechanical_energy(test_args):
    result = mechanical_energy_def.calculate_mechanical_energy(test_args.K, test_args.P)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, mechanical_energy_def.definition_units_SI).evalf(2)
    assert_approx(result_energy, 6)


def test_bad_energy(test_args):
    Eb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        mechanical_energy_def.calculate_mechanical_energy(Eb, test_args.P)
    with raises(TypeError):
        mechanical_energy_def.calculate_mechanical_energy(100, test_args.P)
    with raises(errors.UnitsError):
        mechanical_energy_def.calculate_mechanical_energy(test_args.K, Eb)
    with raises(TypeError):
        mechanical_energy_def.calculate_mechanical_energy(test_args.K, 100)
