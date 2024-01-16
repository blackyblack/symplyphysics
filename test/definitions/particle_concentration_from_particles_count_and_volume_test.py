from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI, dimensionless,
)
from symplyphysics.definitions import particle_concentration_from_particles_count_and_volume as concentration_law


@fixture(name="test_args")
def test_args_fixture():
    N = Quantity(1 * dimensionless)
    V = Quantity(3 * units.meter**3)
    Args = namedtuple("Args", ["N", "V"])
    return Args(N=N, V=V)


def test_basic_concentration(test_args):
    result = concentration_law.calculate_concentration(test_args.N, test_args.V)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, dimensionless / units.volume)
    result_concentration = convert_to(result, concentration_law.definition_units_SI).evalf(2)
    assert result_concentration == approx(0.3333, 0.01)


def test_bad_count_of_particles(test_args):
    Nb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        concentration_law.calculate_concentration(Nb, test_args.V)
    # with raises(TypeError):
    #     concentration_law.calculate_concentration(100, test_args.V)


def test_bad_volume(test_args):
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        concentration_law.calculate_concentration(test_args.N, Vb)
    with raises(TypeError):
        concentration_law.calculate_concentration(test_args.N, 100)
