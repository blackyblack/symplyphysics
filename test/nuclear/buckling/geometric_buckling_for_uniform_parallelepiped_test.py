from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.nuclear.buckling import geometric_buckling_for_uniform_parallelepiped as buckling

@fixture
def test_args():
    # cube reactor 20 * 20 * 20
    par_width = Quantity(20 * units.centimeter)
    par_length = Quantity(20 * units.centimeter)
    par_height = Quantity(20 * units.centimeter)
    Args = namedtuple("Args", ["a", "b", "c"])
    return Args(a=par_width, b=par_length, c=par_height)

def test_basic_geometric_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.a, test_args.b, test_args.c)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**-2)
    result_geometric_buckling = convert_to(result, units.centimeter**-2).subs(units.centimeter, 1).evalf(2)
    assert result_geometric_buckling == approx(0.074, 0.01)

def test_bad_width(test_args):
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(ab, test_args.b, test_args.c)
    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(100, test_args.b, test_args.c)

def test_bad_length(test_args):
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(test_args.a, bb, test_args.c)
    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(test_args.a, 100, test_args.c)

def test_bad_height(test_args):
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(test_args.a, test_args.b, cb)
    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(test_args.a, test_args.b, 100)
