from collections import namedtuple
from pytest import fixture, raises
from sympy import S
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
    dimensionless,
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_infinite_multiplication_factor_diffusion_area as buckling


@fixture(name="test_args")
def test_args_fixture():
    infinite_multiplication_factor = 1.0247
    # critical reactor
    effective_multiplication_factor = 1
    diffusion_area = Quantity(60.117 * units.centimeter**2)
    Args = namedtuple("Args", ["k_inf", "k_eff", "L2"])
    return Args(k_inf=infinite_multiplication_factor,
        k_eff=effective_multiplication_factor,
        L2=diffusion_area)


def test_basic_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.k_inf, test_args.k_eff,
        test_args.L2)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.area)
    result_buckling = convert_to(result, 1 / units.centimeter**2).evalf(4)
    assert_approx(result_buckling, 0.0004109)


def test_zero_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(1, 1, test_args.L2)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, dimensionless)
    result_buckling = convert_to(result, S.One).evalf(4)
    assert_approx(result_buckling, 0)


def test_bad_diffusion_area(test_args):
    Lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(test_args.k_inf, test_args.k_eff, Lb)
    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(test_args.k_inf, test_args.k_eff, 100)
