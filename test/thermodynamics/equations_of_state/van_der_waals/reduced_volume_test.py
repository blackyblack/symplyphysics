from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import reduced_volume

# Description
## The volume of a van der Waals gas is 1 m**3, the value of critical volume is 0.5 m**3. Then the reduced
## volume of the gas is 2.

Args = namedtuple("Args", "v vc")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v = Quantity(1 * units.meter**3)
    vc = Quantity(0.5 * units.meter**3)
    return Args(v=v, vc=vc)


def test_law(test_args: Args) -> None:
    result = reduced_volume.calculate_reduced_volume(test_args.v, test_args.vc)
    assert_equal(result, 2)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        reduced_volume.calculate_reduced_volume(vb, test_args.vc)
    with raises(TypeError):
        reduced_volume.calculate_reduced_volume(100, test_args.vc)
    with raises(errors.UnitsError):
        reduced_volume.calculate_reduced_volume(test_args.v, vb)
    with raises(TypeError):
        reduced_volume.calculate_reduced_volume(test_args.v, 100)
