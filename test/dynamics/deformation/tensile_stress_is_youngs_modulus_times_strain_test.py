from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import (
    tensile_stress_is_youngs_modulus_times_strain as stress_strain_law
)

# Description
## An object is being compressed, its initial length is 1 m, and it became 2 mm shorter.
## The object's Young's modulus is 10 Pa. The stress on the object amounts to 2e-2 Pa.


@fixture(name="test_args")
def test_args_fixture():
    L = Quantity(1.0 * units.meter)
    dL = Quantity(-2.0 * units.millimeter)
    E = Quantity(10.0 * units.pascal)
    Args = namedtuple("Args", "L dL E")
    return Args(L=L, dL=dL, E=E)


def test_law(test_args):
    result = stress_strain_law.calculate_tensile_stress(test_args.E, test_args.dL, test_args.L)
    assert_equal(result, 2e-2 * units.pascal)


def test_bad_modulus(test_args):
    Eb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        stress_strain_law.calculate_tensile_stress(Eb, test_args.dL, test_args.L)
    with raises(TypeError):
        stress_strain_law.calculate_tensile_stress(100, test_args.dL, test_args.L)


def test_bad_length(test_args):
    Lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        stress_strain_law.calculate_tensile_stress(test_args.E, Lb, test_args.L)
    with raises(TypeError):
        stress_strain_law.calculate_tensile_stress(test_args.E, 100, test_args.L)
    with raises(errors.UnitsError):
        stress_strain_law.calculate_tensile_stress(test_args.E, test_args.dL, Lb)
    with raises(TypeError):
        stress_strain_law.calculate_tensile_stress(test_args.E, test_args.dL, 100)
