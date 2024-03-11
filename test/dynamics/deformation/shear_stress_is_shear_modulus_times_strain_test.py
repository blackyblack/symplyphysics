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
## Block of wood (G = 4 GPa) is undergoing shear stress. The transverse displacement of the
## wood's surface is 0.5 mm. The length of the face of the block which the shear force is applied
## to is 40 cm. The shear stress amounts to 5 MPa.

Args = namedtuple("Args", "g dx l")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    g = Quantity(4 * prefixes.giga * units.pascal)
    dx = Quantity(0.5 * units.millimeter)
    l = Quantity(40 * units.centimeter)
    return Args(g=g, dx=dx, l=l)


def test_law(test_args: Args) -> None:
    result = shear_stress_law.calculate_shear_stress(test_args.g, test_args.dx, test_args.l)
    assert_equal(result, 5 * prefixes.mega * units.pascal)


def test_bad_shear_modulus(test_args: Args) -> None:
    gb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        shear_stress_law.calculate_shear_stress(gb, test_args.dx, test_args.l)
    with raises(TypeError):
        shear_stress_law.calculate_shear_stress(100, test_args.dx, test_args.l)


def test_bad_length(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        shear_stress_law.calculate_shear_stress(test_args.g, lb, test_args.l)
    with raises(TypeError):
        shear_stress_law.calculate_shear_stress(test_args.g, 100, test_args.l)
    with raises(errors.UnitsError):
        shear_stress_law.calculate_shear_stress(test_args.g, test_args.dx, lb)
    with raises(TypeError):
        shear_stress_law.calculate_shear_stress(test_args.g, test_args.dx, 100)
