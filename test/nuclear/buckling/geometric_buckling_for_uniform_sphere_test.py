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
from symplyphysics.laws.nuclear.buckling import geometric_buckling_for_uniform_sphere as buckling


@fixture(name="test_args")
def test_args_fixture():
    sphere_radius = Quantity(20 * units.centimeter)
    Args = namedtuple("Args", ["R"])
    return Args(R=sphere_radius)


def test_basic_geometric_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.R)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.area)
    result_geometric_buckling = convert_to(result, 1 / units.centimeter**2).evalf(4)
    assert_approx(result_geometric_buckling, 0.02467)


def test_bad_sphere_radius():
    Rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(Rb)
    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(100)
