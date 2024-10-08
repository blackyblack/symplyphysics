from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.optics import optical_path_length_from_geometrical_path_length_and_refractive_index as optical_path

Args = namedtuple("Args", ["geometric_path", "refractive_index"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    geometric_path = Quantity(0.001 * units.meter)
    refractive_index = 1.5
    return Args(
        geometric_path=geometric_path,
        refractive_index=refractive_index,
    )


def test_basic_optical_path(test_args: Args) -> None:
    result = optical_path.calculate_optical_path(test_args.geometric_path,
        test_args.refractive_index)
    assert_equal(result, 0.0015 * units.meter)


def test_bad_geometric_path(test_args: Args) -> None:
    bad_geometric_path = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        optical_path.calculate_optical_path(bad_geometric_path, test_args.refractive_index)
    with raises(TypeError):
        optical_path.calculate_optical_path(100, test_args.refractive_index)


def test_bad_refractive_index(test_args: Args) -> None:
    bad_refractive_index = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        optical_path.calculate_optical_path(test_args.geometric_path, bad_refractive_index)
