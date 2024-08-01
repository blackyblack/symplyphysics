from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import (tensile_stress_is_youngs_modulus_times_strain
    as stress_strain_law)

# Description
## An object is being compressed, its strain amounts to -0.002. The object's Young's
## modulus is 10 Pa. The stress on the object amounts to -0.02 Pa.

Args = namedtuple("Args", "e s")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = Quantity(10.0 * units.pascal)
    s = -0.002
    return Args(e=e, s=s)


def test_law(test_args: Args) -> None:
    result = stress_strain_law.calculate_tensile_stress(test_args.e, test_args.s)
    assert_equal(result, -0.02 * units.pascal)


def test_bad_modulus(test_args: Args) -> None:
    Eb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        stress_strain_law.calculate_tensile_stress(Eb, test_args.s)
    with raises(TypeError):
        stress_strain_law.calculate_tensile_stress(100, test_args.s)


def test_bad_strain(test_args: Args) -> None:
    sb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        stress_strain_law.calculate_tensile_stress(test_args.e, sb)
