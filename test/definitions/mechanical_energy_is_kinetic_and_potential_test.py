from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import mechanical_energy_is_kinetic_and_potential as mechanical_energy_def


@fixture
def test_args():
    K = Quantity(1 * units.joule)
    P = Quantity(5 * units.joule)
    Args = namedtuple("Args", ["K", "P"])
    return Args(K=K, P=P)


def test_basic_mechanical_energy(test_args):
    result = mechanical_energy_def.calculate_mechanical_energy(test_args.K, test_args.P)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_momentum = convert_to(result,
        mechanical_energy_def.definition_units_SI).subs(units.joule, 1).evalf(2)
    assert result_momentum == approx(6.0, 0.001)


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
