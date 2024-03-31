from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import wave_equation_general_solution_in_one_dimension as wave_law

# Description
## A wave is traveling in a stretched string with maximum amplitude of 5 cm. When wave phase
## is 2.3 rad, the displacement of the string is -3.33 cm.

Args = namedtuple("Args", "a psi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(5.0 * units.centimeter)
    psi = 2.3
    return Args(a=a, psi=psi)


def test_law(test_args: Args) -> None:
    result = wave_law.calculate_displacement(test_args.a, test_args.psi)
    assert_equal(result, -3.33 * units.centimeter)


def test_bad_wave_phase(test_args: Args) -> None:
    psi_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_law.calculate_displacement(test_args.a, psi_bad)
