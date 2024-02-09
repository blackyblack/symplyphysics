from collections import namedtuple
from pytest import fixture, raises
from sympy import cos, pi
from sympy.vector import CoordSys3D
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear import diffusion_equation_from_neutron_flux as diffusion_equation

Args = namedtuple("Args", ["f", "v", "Sf", "Sa", "D"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # cube reactor with side = 1 meter
    cartesian_coordinates = CoordSys3D("cartesian_coordinates")
    # Makes linter happy
    x = getattr(cartesian_coordinates, "x")
    y = getattr(cartesian_coordinates, "y")
    z = getattr(cartesian_coordinates, "z")
    cube_side = Quantity(1 * units.meter)
    unit_length = Quantity(1 * units.meter)
    neutron_flux_unit = Quantity(1 / units.meter**2 / units.second)
    neutron_flux = neutron_flux_unit * (cos(pi / cube_side * x * unit_length) *
        cos(pi / cube_side * y * unit_length) * cos(pi / cube_side * z * unit_length))
    neutrons_per_fission = 1
    macro_fission_cross_section = Quantity(0.006 / units.centimeter)
    macro_abs_cross_section = Quantity(0.0025 / units.centimeter)
    diffusion_coefficient = Quantity(2 * units.centimeter)
    return Args(f=neutron_flux,
        v=neutrons_per_fission,
        Sf=macro_fission_cross_section,
        Sa=macro_abs_cross_section,
        D=diffusion_coefficient)


def test_basic_multiplication_factor(test_args: Args) -> None:
    result = diffusion_equation.calculate_multiplication_factor(test_args.f, test_args.v,
        test_args.Sf, test_args.Sa, test_args.D)
    assert_equal(result, 0.712)


def test_bad_diffusion_coefficient(test_args: Args) -> None:
    Db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_equation.calculate_multiplication_factor(test_args.f, test_args.v, test_args.Sf,
            test_args.Sa, Db)
    with raises(TypeError):
        diffusion_equation.calculate_multiplication_factor(test_args.f, test_args.v, test_args.Sf,
            test_args.Sa, 100)


def test_bad_macroscopic_cross_section(test_args: Args) -> None:
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_equation.calculate_multiplication_factor(test_args.f, test_args.v, Sb,
            test_args.Sa, test_args.D)
    with raises(TypeError):
        diffusion_equation.calculate_multiplication_factor(test_args.f, test_args.v, 100,
            test_args.Sa, test_args.D)
    with raises(errors.UnitsError):
        diffusion_equation.calculate_multiplication_factor(test_args.f, test_args.v, test_args.Sf,
            Sb, test_args.D)
    with raises(TypeError):
        diffusion_equation.calculate_multiplication_factor(test_args.f, test_args.v, test_args.Sf,
            100, test_args.D)
