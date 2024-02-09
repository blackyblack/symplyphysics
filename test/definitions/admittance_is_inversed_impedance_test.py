from collections import namedtuple
from sympy import I
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import admittance_is_inversed_impedance as admittance_def

# Description
## If the dipole has 2i Ohm impedance, it should have -0.5i Siemens admittance. No external calculators were used for such computation.

Args = namedtuple("Args", ["Z"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    Z = Quantity(2 * I * units.ohm)
    return Args(Z=Z)


def test_basic_admittance(test_args: Args) -> None:
    result = admittance_def.calculate_admittance(test_args.Z)
    assert_equal(result, -0.5 * I * units.siemens)


def test_bad_impedance() -> None:
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        admittance_def.calculate_admittance(Rb)
    with raises(TypeError):
        admittance_def.calculate_admittance(100)
