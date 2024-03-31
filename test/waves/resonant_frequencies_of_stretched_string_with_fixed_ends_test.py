from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import (
    resonant_frequencies_of_stretched_string_with_fixed_ends as frequencies_law,)

# Description
## The resonant frequency of the second harmonic of the wave traveling with a phase speed of 10 m/s
## along a stretched string of a length of 40 cm is 25 Hz.

Args = namedtuple("Args", "n v l")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n = 2
    v = Quantity(10 * units.meter / units.second)
    l = Quantity(40 * units.centimeter)
    return Args(n=n, v=v, l=l)


def test_law(test_args: Args) -> None:
    result = frequencies_law.calculate_resonant_frequency(test_args.n, test_args.v, test_args.l)
    assert_equal(result, 25 * units.hertz)


def test_bad_harmonic_number(test_args: Args) -> None:
    nb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        frequencies_law.calculate_resonant_frequency(nb, test_args.v, test_args.l)


def test_bad_phase_velocity(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        frequencies_law.calculate_resonant_frequency(test_args.n, vb, test_args.l)
    with raises(TypeError):
        frequencies_law.calculate_resonant_frequency(test_args.n, 100, test_args.l)


def test_bad_string_length(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        frequencies_law.calculate_resonant_frequency(test_args.n, test_args.v, lb)
    with raises(TypeError):
        frequencies_law.calculate_resonant_frequency(test_args.n, test_args.v, 100)
