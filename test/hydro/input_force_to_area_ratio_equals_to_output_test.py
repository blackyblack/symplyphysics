from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)

from symplyphysics.laws.hydro import input_force_to_area_ratio_equals_to_output as hydraulic


@fixture(name="test_args")
def test_args_fixture():
    first_force = Quantity(2000 * units.newton)
    first_area = Quantity(10 * units.meter**2)
    second_area = Quantity(0.1 * units.meter**2)
    Args = namedtuple("Args", ["first_force", "first_area", "second_area"])
    return Args(first_force=first_force, first_area=first_area, second_area=second_area)


def test_basic_force(test_args):
    result = hydraulic.calculate_output_force(test_args.first_force, test_args.first_area,
        test_args.second_area)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force = convert_to(result, units.newton).evalf(5)
    assert result_force == approx(20, 0.001)


def test_bad_area(test_args):
    ad = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hydraulic.calculate_output_force(test_args.first_force, ad, test_args.second_area)
    with raises(errors.UnitsError):
        hydraulic.calculate_output_force(test_args.first_force, test_args.first_area, ad)
    with raises(TypeError):
        hydraulic.calculate_output_force(test_args.first_force, 100, test_args.second_area)
    with raises(TypeError):
        hydraulic.calculate_output_force(test_args.first_force, test_args.first_area, 100)


def test_bad_force(test_args):
    ag = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hydraulic.calculate_output_force(ag, test_args.first_area, test_args.second_area)
    with raises(TypeError):
        hydraulic.calculate_output_force(100, test_args.first_area, test_args.second_area)
