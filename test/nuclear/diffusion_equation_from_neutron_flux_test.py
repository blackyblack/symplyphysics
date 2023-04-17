from collections import namedtuple
from pytest import approx, fixture, raises

from sympy.vector import CoordSys3D
from symplyphysics import (
    units, errors, cos, pi
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.nuclear import diffusion_equation_from_neutron_flux as diffusion_equation

@fixture
def test_args():
    # cube reactor with side = 1 meter
    cartesian_coordinates = CoordSys3D("cartesian_coordinates")
    cube_side = Quantity(units.length, 1 * units.meter)
    neutron_flux = (cos((pi / cube_side) * cartesian_coordinates.x) *
        cos((pi / cube_side) * cartesian_coordinates.y) *
        cos((pi / cube_side) * cartesian_coordinates.z))
    neutrons_per_fission = 1
    macro_fission_cross_section = Quantity(1 / units.length, 0.006 / units.centimeter)
    macro_abs_cross_section = Quantity(1 / units.length, 0.0025 / units.centimeter)
    diffusion_coefficient = Quantity(units.length, 2 * units.centimeter)
    Args = namedtuple("Args", ["f", "v", "Sf", "Sa", "D"])
    return Args(f=neutron_flux, v=neutrons_per_fission,
        Sf=macro_fission_cross_section, Sa=macro_abs_cross_section, D=diffusion_coefficient)

def test_basic_multiplication_factor(test_args):
    result = diffusion_equation.calculate_multiplication_factor(
        test_args.f, test_args.v, test_args.Sf, test_args.Sa, test_args.D)
    assert result == approx(0.712, 0.01)

def test_bad_diffusion_coefficient(test_args):
    Db = Quantity(units.time)
    with raises(errors.UnitsError):
        diffusion_equation.calculate_multiplication_factor(
            test_args.f, test_args.v, test_args.Sf, test_args.Sa, Db)
    with raises(TypeError):
        diffusion_equation.calculate_multiplication_factor(
            test_args.f, test_args.v, test_args.Sf, test_args.Sa, 100)

def test_bad_macroscopic_cross_section(test_args):
    Sb = Quantity(units.time)
    with raises(errors.UnitsError):
        diffusion_equation.calculate_multiplication_factor(
            test_args.f, test_args.v, Sb, test_args.Sa, test_args.D)
    with raises(TypeError):
        diffusion_equation.calculate_multiplication_factor(
            test_args.f, test_args.v, 100, test_args.Sa, test_args.D)
    with raises(errors.UnitsError):
        diffusion_equation.calculate_multiplication_factor(
            test_args.f, test_args.v, test_args.Sf, Sb, test_args.D)
    with raises(TypeError):
        diffusion_equation.calculate_multiplication_factor(
            test_args.f, test_args.v, test_args.Sf, 100, test_args.D)
