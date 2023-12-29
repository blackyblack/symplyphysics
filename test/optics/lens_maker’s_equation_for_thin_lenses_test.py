from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)

from symplyphysics.laws.optics import lens_makers_equation_for_thin_lenses as makers_equation


@fixture(name="test_args")
def test_args_fixture():
    n = 1.52
    radius_1 = Quantity(0.1 * units.meter)
    radius_2 = Quantity(0.2 * units.meter)
    Args = namedtuple("Args", ["n", "radius_1", "radius_2"])
    return Args(n=n,
        radius_1=radius_1,
        radius_2=radius_2,
               )


def test_basic_length(test_args):
    result = makers_equation.calculate_length(test_args.n, test_args.radius_1, test_args.radius_2)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_length = convert_to(result, units.meter).evalf(5)
    assert result_length == approx(0.38493, 0.00001)


def test_bad_radius(test_args):
    ag = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        makers_equation.calculate_length(ag, test_args.radius_1, test_args.radius_2)
    with raises(errors.UnitsError):
        makers_equation.calculate_length(test_args.n, ag, test_args.radius_2)
    with raises(errors.UnitsError):
        makers_equation.calculate_length(test_args.n, test_args.radius_1, ag)
    with raises(TypeError):
        makers_equation.calculate_length(test_args.n, 100, test_args.radius_2)
    with raises(TypeError):
        makers_equation.calculate_length(test_args.n, test_args.radius_1, 100)

