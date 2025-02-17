from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import volumetric_coefficient_of_thermal_expansion as coefficient_def

# Description
## A body is heated up and changes its volume from 1.00 m**3 at temperature T = 400 K to 1.001 m**3 at
## T = 401 K. Then its volumetric coefficient of thermal expansion is 9.9e-3 1/K (at T = 401 K).

Args = namedtuple("Args", "v0 v1 t0 t1")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v0 = Quantity(1.00 * units.meter**3)
    v1 = Quantity(1.01 * units.meter**3)
    t0 = Quantity(400 * units.kelvin)
    t1 = Quantity(401 * units.kelvin)
    return Args(v0=v0, v1=v1, t0=t0, t1=t1)


def test_law(test_args: Args) -> None:
    result = coefficient_def.calculate_volumetric_expansion_coefficient(
        test_args.v0, test_args.v1, test_args.t0, test_args.t1)
    assert_equal(result, 9.9e-3 / units.kelvin, relative_tolerance=6e-3)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_def.calculate_volumetric_expansion_coefficient(vb, test_args.v1, test_args.t0,
            test_args.t1)
    with raises(TypeError):
        coefficient_def.calculate_volumetric_expansion_coefficient(100, test_args.v1, test_args.t0,
            test_args.t1)
    with raises(errors.UnitsError):
        coefficient_def.calculate_volumetric_expansion_coefficient(test_args.v0, vb, test_args.t0,
            test_args.t1)
    with raises(TypeError):
        coefficient_def.calculate_volumetric_expansion_coefficient(test_args.v0, 100, test_args.t0,
            test_args.t1)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_def.calculate_volumetric_expansion_coefficient(test_args.v0, test_args.v1, tb,
            test_args.t1)
    with raises(TypeError):
        coefficient_def.calculate_volumetric_expansion_coefficient(test_args.v0, test_args.v1, 100,
            test_args.t1)
    with raises(errors.UnitsError):
        coefficient_def.calculate_volumetric_expansion_coefficient(test_args.v0, test_args.v1,
            test_args.t0, tb)
    with raises(TypeError):
        coefficient_def.calculate_volumetric_expansion_coefficient(test_args.v0, test_args.v1,
            test_args.t0, 100)
