from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, units, Quantity, errors
from symplyphysics.laws.thermodynamics import (
    isochoric_and_isobaric_heat_capacities_of_homogeneous_substance as mayers_relation,)

# Description
## For a homogeneous substance with isobaric heat capacity C_p = 5 J/K, expansion coefficient
## alpha = 2e-4 1/K, isothermal compressibility beta_T = 4.5e-10 1/Pa, of volume V = 0.1 liter
## and at T = 300 K, the isochoric heat capacity amounts to approximately 2.33 J/K.

Args = namedtuple("Args", "cp v t a b")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cp = Quantity(5 * units.joule / units.kelvin)
    v = Quantity(0.1 * units.liter)
    t = Quantity(300 * units.kelvin)
    a = Quantity(2e-4 / units.kelvin)
    b = Quantity(4.5e-10 / units.pascal)
    return Args(cp=cp, v=v, t=t, a=a, b=b)


def test_law(test_args: Args) -> None:
    result = mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, test_args.v,
        test_args.t, test_args.a, test_args.b)
    assert_equal(result, 2.33 * units.joule / units.kelvin, relative_tolerance=2e-3)


def test_bad_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        mayers_relation.calculate_isochoric_heat_capacity(cb, test_args.v, test_args.t, test_args.a,
            test_args.b)
    with raises(TypeError):
        mayers_relation.calculate_isochoric_heat_capacity(100, test_args.v, test_args.t,
            test_args.a, test_args.b)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, vb, test_args.t,
            test_args.a, test_args.b)
    with raises(TypeError):
        mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, 100, test_args.t,
            test_args.a, test_args.b)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, test_args.v, tb,
            test_args.a, test_args.b)
    with raises(TypeError):
        mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, test_args.v, 100,
            test_args.a, test_args.b)


def test_bad_expansion_coefficient(test_args: Args) -> None:
    ab = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, test_args.v, test_args.t,
            ab, test_args.b)
    with raises(TypeError):
        mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, test_args.v, test_args.t,
            100, test_args.b)


def test_bad_compressibility(test_args: Args) -> None:
    bb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, test_args.v, test_args.t,
            test_args.a, bb)
    with raises(TypeError):
        mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, test_args.v, test_args.t,
            test_args.a, 100)
