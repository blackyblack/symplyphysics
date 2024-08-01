from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import (
    engineering_normal_strain_is_total_deformation_over_initial_dimension as law,
)

Args = namedtuple("Args", "dl l")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dl = Quantity(3 * units.centimeter)
    l = Quantity(200 * units.meter)
    return Args(dl=dl, l=l)


def test_law(test_args: Args) -> None:
    result = law.calculate_engineering_normal_strain(test_args.dl, test_args.l)
    assert_equal(result, 1.5e-4)


def test_bad_length(test_args: Args) -> None:
    lb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_engineering_normal_strain(lb, test_args.l)
    with raises(TypeError):
        law.calculate_engineering_normal_strain(100, test_args.l)
    with raises(errors.UnitsError):
        law.calculate_engineering_normal_strain(test_args.dl, lb)
    with raises(TypeError):
        law.calculate_engineering_normal_strain(test_args.dl, 100)
