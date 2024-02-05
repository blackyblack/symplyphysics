from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.nuclear.buckling import material_buckling_from_macroscopic_fission_cross_section_diffusion_coefficient as material_buckling


@fixture(name="test_args")
def test_args_fixture():
    neutrons_per_fission = 2.6
    macro_fission_cross_section = Quantity(1.482 / units.centimeter)
    macro_abs_cross_section = Quantity(3.108 / units.centimeter)
    diffusion_coefficient = Quantity(31.782 * units.meter)
    Args = namedtuple("Args", ["v", "Sf", "Sa", "D"])
    return Args(v=neutrons_per_fission,
        Sf=macro_fission_cross_section,
        Sa=macro_abs_cross_section,
        D=diffusion_coefficient)


def test_basic_buckling(test_args):
    result = material_buckling.calculate_buckling(test_args.v, test_args.Sf, test_args.Sa,
        test_args.D)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.area)
    result_buckling = convert_to(result, 1 / units.meter**2).evalf(4)
    assert_approx(result_buckling, 2.345)


def test_bad_macroscopic_cross_section(test_args):
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        material_buckling.calculate_buckling(test_args.v, Sb, test_args.Sa, test_args.D)
    with raises(TypeError):
        material_buckling.calculate_buckling(test_args.v, 100, test_args.Sa, test_args.D)
    with raises(errors.UnitsError):
        material_buckling.calculate_buckling(test_args.v, test_args.Sf, Sb, test_args.D)
    with raises(TypeError):
        material_buckling.calculate_buckling(test_args.v, test_args.Sf, 100, test_args.D)


def test_bad_diffusion_coefficient(test_args):
    Db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        material_buckling.calculate_buckling(test_args.v, test_args.Sf, test_args.Sa, Db)
    with raises(TypeError):
        material_buckling.calculate_buckling(test_args.v, test_args.Sf, test_args.Sa, 100)
