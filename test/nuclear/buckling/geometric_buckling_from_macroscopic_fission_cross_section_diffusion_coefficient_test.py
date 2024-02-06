from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_macroscopic_fission_cross_section_diffusion_coefficient as buckling


@fixture(name="test_args")
def test_args_fixture():
    neutrons_per_fission = 2.6
    # critical reactor
    effective_multiplication_factor = 1
    macro_fission_cross_section = Quantity(1.482 / units.centimeter)
    macro_abs_cross_section = Quantity(3.108 / units.centimeter)
    diffusion_coefficient = Quantity(31.782 * units.meter)
    Args = namedtuple("Args", ["v", "k", "Sf", "Sa", "D"])
    return Args(v=neutrons_per_fission,
        k=effective_multiplication_factor,
        Sf=macro_fission_cross_section,
        Sa=macro_abs_cross_section,
        D=diffusion_coefficient)


def test_basic_buckling(test_args):
    result = buckling.calculate_buckling(test_args.v, test_args.k, test_args.Sf, test_args.Sa,
        test_args.D)
    assert_equal(result, 2.345 / units.meter**2)


def test_bad_macroscopic_cross_section(test_args):
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_buckling(test_args.v, test_args.k, Sb, test_args.Sa, test_args.D)
    with raises(TypeError):
        buckling.calculate_buckling(test_args.v, test_args.k, 100, test_args.Sa, test_args.D)
    with raises(errors.UnitsError):
        buckling.calculate_buckling(test_args.v, test_args.k, test_args.Sf, Sb, test_args.D)
    with raises(TypeError):
        buckling.calculate_buckling(test_args.v, test_args.k, test_args.Sf, 100, test_args.D)


def test_bad_diffusion_coefficient(test_args):
    Db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_buckling(test_args.v, test_args.k, test_args.Sf, test_args.Sa, Db)
    with raises(TypeError):
        buckling.calculate_buckling(test_args.v, test_args.k, test_args.Sf, test_args.Sa, 100)
