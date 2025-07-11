from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.hydro import force_to_area_ratio_in_hydraulic_press as hydraulic

Args = namedtuple("Args", ["first_force", "first_area", "second_area"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    first_force = Quantity(2000 * units.newton)
    first_area = Quantity(10 * units.meter**2)
    second_area = Quantity(0.1 * units.meter**2)
    return Args(first_force=first_force, first_area=first_area, second_area=second_area)


def test_basic_force(test_args: Args) -> None:
    result = hydraulic.calculate_output_force(test_args.first_force, test_args.first_area,
        test_args.second_area)
    assert_equal(result, 20 * units.newton)


def test_bad_area(test_args: Args) -> None:
    ad = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hydraulic.calculate_output_force(test_args.first_force, ad, test_args.second_area)
    with raises(errors.UnitsError):
        hydraulic.calculate_output_force(test_args.first_force, test_args.first_area, ad)
    with raises(TypeError):
        hydraulic.calculate_output_force(test_args.first_force, 100, test_args.second_area)
    with raises(TypeError):
        hydraulic.calculate_output_force(test_args.first_force, test_args.first_area, 100)


def test_bad_force(test_args: Args) -> None:
    ag = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hydraulic.calculate_output_force(ag, test_args.first_area, test_args.second_area)
    with raises(TypeError):
        hydraulic.calculate_output_force(100, test_args.first_area, test_args.second_area)
