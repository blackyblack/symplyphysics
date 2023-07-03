from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import volume_number_density


@fixture(name="test_args")
def test_args_fixture_fixture():
    objects = 100
    V = Quantity(1 * units.meter**3)
    Args = namedtuple("Args", ["o", "V"])
    return Args(o=objects, V=V)


def test_basic_density(test_args):
    result = volume_number_density.calculate_number_density(test_args.o, test_args.V)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.volume)
    result_density = convert_to(result, volume_number_density.definition_units_SI).evalf(2)
    assert result_density == approx(100, 0.01)


def test_bad_volume(test_args):
    Vb = Quantity(1 * units.length)
    with raises(errors.UnitsError):
        volume_number_density.calculate_number_density(test_args.o, Vb)
    with raises(TypeError):
        volume_number_density.calculate_number_density(test_args.o, 100)
