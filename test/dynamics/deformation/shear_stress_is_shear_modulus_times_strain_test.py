from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    prefixes,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import (
    shear_stress_is_shear_modulus_times_strain as shear_stress_law,
)

# Description
## Block of wood (G = 4 GPa) is undergoing shear stress. The shear strain (angle of deformation) is 1.25e-3 rad.
## The shear stress in the block amounts to 5 MPa.

Args = namedtuple("Args", "g gamma")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    g = Quantity(4 * prefixes.giga * units.pascal)
    gamma = Quantity(1.25e-3 * units.radian)
    return Args(g=g, gamma=gamma)


def test_law(test_args: Args) -> None:
    result = shear_stress_law.calculate_shear_stress(test_args.g, test_args.gamma)
    assert_equal(result, 5 * prefixes.mega * units.pascal)


def test_bad_shear_modulus(test_args: Args) -> None:
    gb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        shear_stress_law.calculate_shear_stress(gb, test_args.gamma)
    with raises(TypeError):
        shear_stress_law.calculate_shear_stress(100, test_args.gamma)


def test_bad_angle(test_args: Args) -> None:
    gamma_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        shear_stress_law.calculate_shear_stress(test_args.g, gamma_bad)
