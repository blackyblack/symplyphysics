from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_for_uniform_sphere as buckling


@fixture
def test_args():
    sphere_radius = Quantity(20 * units.centimeter)
    Args = namedtuple("Args", ["R"])
    return Args(R=sphere_radius)


def test_basic_geometric_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.R)
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
                                                     units.length**-2)
    result_geometric_buckling = convert_to(result, units.centimeter**-2).subs(
        units.centimeter, 1).evalf(2)
    assert result_geometric_buckling == approx(0.0246, 0.01)


def test_bad_sphere_radius():
    Rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(Rb)
    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(100)
