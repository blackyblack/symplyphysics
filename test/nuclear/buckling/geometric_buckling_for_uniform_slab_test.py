from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_for_uniform_slab as buckling


@fixture(name="test_args")
def test_args_fixture():
    slab_width = Quantity(200 * units.centimeter)
    Args = namedtuple("Args", ["A"])
    return Args(A=slab_width)


def test_basic_geometric_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.A)
    assert_equal(result, 0.0002467 / units.centimeter**2)


def test_bad_slab_width():
    Ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(Ab)
    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(100)
