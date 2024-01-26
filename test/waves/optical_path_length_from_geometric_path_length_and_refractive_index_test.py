from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
    dimensionless
)

from symplyphysics.laws.waves import optical_path_length_from_geometric_path_length_and_refractive_index as optical_path


@fixture(name="test_args")
def test_args_fixture():
    geometric_path = Quantity(0.001 * units.meter)
    refractive_index = 1.5
    Args = namedtuple("Args", ["geometric_path", "refractive_index"])
    return Args(
        geometric_path=geometric_path,
        refractive_index=refractive_index,
               )


def test_basic_optical_path(test_args):
    result = optical_path.calculate_optical_path(test_args.geometric_path, test_args.refractive_index)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_optical_path = convert_to(result, units.meter.evalf(5))
    assert result_optical_path == approx(0.0015, 0.0001)


def test_bad_geometric_path(test_args):
    bad_geometric_path = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        optical_path.calculate_optical_path(bad_geometric_path, test_args.refractive_index)
    with raises(TypeError):
        optical_path.calculate_optical_path(100, test_args.refractive_index)


def test_bad_refractive_index(test_args):
    bad_refractive_index = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        optical_path.calculate_optical_path(test_args.geometric_path, bad_refractive_index)
