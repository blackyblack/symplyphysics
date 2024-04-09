from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.hydro import hydraulic_stress_is_bulk_modulus_times_strain as stress_law

# Description
## A steel object (B = 160 GPa) is undergoing bulk compression in a fluid.
## Its initial volume was 20 cm**3 and it shrunk by 0.5 mm**3. Then the hydraulic
## pressure on the object amounts to 4 MPa.

Args = namedtuple("Args", "b dv v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    b = Quantity(160 * prefixes.giga * units.pascal)
    dv = Quantity(0.5 * units.millimeter**3)
    v = Quantity(20 * units.centimeter**3)
    return Args(b=b, dv=dv, v=v)


def test_law(test_args: Args) -> None:
    result = stress_law.calculate_hydraulic_stress(test_args.b, test_args.dv, test_args.v)
    assert_equal(result, 4 * prefixes.mega * units.pascal)


def test_bad_bulk_modulus(test_args: Args) -> None:
    bb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        stress_law.calculate_hydraulic_stress(bb, test_args.dv, test_args.v)
    with raises(TypeError):
        stress_law.calculate_hydraulic_stress(100, test_args.dv, test_args.v)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        stress_law.calculate_hydraulic_stress(test_args.b, vb, test_args.v)
    with raises(TypeError):
        stress_law.calculate_hydraulic_stress(test_args.b, 100, test_args.v)
    with raises(errors.UnitsError):
        stress_law.calculate_hydraulic_stress(test_args.b, test_args.dv, vb)
    with raises(TypeError):
        stress_law.calculate_hydraulic_stress(test_args.b, test_args.dv, 100)
