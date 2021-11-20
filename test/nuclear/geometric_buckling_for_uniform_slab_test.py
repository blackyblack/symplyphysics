from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.nuclear import geometric_buckling_for_uniform_slab as buckling

@fixture
def test_args():
    slab_width = units.Quantity('slab_width')
    SI.set_quantity_dimension(slab_width, units.length)
    SI.set_quantity_scale_factor(slab_width, 200 * units.centimeter)

    Args = namedtuple('Args', ['A'])
    return Args(A = slab_width)

def test_basic_geometric_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.A)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**-2)

    result_geometric_buckling = convert_to(result, units.centimeter**-2).subs(units.centimeter, 1).evalf(2)
    assert result_geometric_buckling == approx(0.000246, 0.01)

def test_bad_slab_width():
    Ab = units.Quantity('Ab')
    SI.set_quantity_dimension(Ab, units.time)
    SI.set_quantity_scale_factor(Ab, 3 * units.second)

    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(Ab)

    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(100)
