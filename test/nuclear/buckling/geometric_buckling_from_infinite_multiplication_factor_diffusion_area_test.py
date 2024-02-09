from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_infinite_multiplication_factor_diffusion_area as buckling

Args = namedtuple("Args", ["k_inf", "k_eff", "L2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    infinite_multiplication_factor = 1.0247
    # critical reactor
    effective_multiplication_factor = 1
    diffusion_area = Quantity(60.117 * units.centimeter**2)
    return Args(k_inf=infinite_multiplication_factor,
        k_eff=effective_multiplication_factor,
        L2=diffusion_area)


def test_basic_buckling(test_args: Args) -> None:
    result = buckling.calculate_geometric_buckling_squared(test_args.k_inf, test_args.k_eff,
        test_args.L2)
    assert_equal(result, 0.0004109 / units.centimeter**2)


def test_zero_buckling(test_args: Args) -> None:
    result = buckling.calculate_geometric_buckling_squared(1, 1, test_args.L2)
    assert_equal(result, 0)


def test_bad_diffusion_area(test_args: Args) -> None:
    Lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(test_args.k_inf, test_args.k_eff, Lb)
    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(test_args.k_inf, test_args.k_eff, 100)
