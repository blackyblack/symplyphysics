from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.equations_of_state import (
    dimensionless_van_den_waals_equation as state_equation_law,
)

# Description
## The value of reduced pressure for any van der Waals fluid is p* = 0.176 at reduced
## volume V* = 3.5 and reduced temperature T* = 0.50

Args = namedtuple("Args", "vr tr")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    vr = 3.5
    tr = 0.5
    return Args(vr=vr, tr=tr)


def test_law(test_args: Args) -> None:
    result = state_equation_law.calculate_reduced_pressure(test_args.vr, test_args.tr)
    assert_equal(result, 0.176)


def test_bad_reduced_quantity(test_args: Args) -> None:
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        state_equation_law.calculate_reduced_pressure(rb, test_args.tr)    
    with raises(errors.UnitsError):
        state_equation_law.calculate_reduced_pressure(test_args.vr, rb)
