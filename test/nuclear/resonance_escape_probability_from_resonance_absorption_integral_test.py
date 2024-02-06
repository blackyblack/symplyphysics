from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear import resonance_escape_probability_from_resonance_absorption_integral as resonance_escape


@fixture(name="test_args")
def test_args_fixture():
    # U-238 atomic number density = 0.984 atoms per (barn * cm)
    atomic_number_density_abs = Quantity(0.984e24 / units.centimeter**3)
    # UO2 fuel rod with diameter = 1 cm has effective resonance integral = 20.5 barns
    resonance_integral = Quantity(2.05e-23 * units.centimeter**2)
    # U-238 with carbon moderator average lethargy decrement = 0.1573
    average_lethargy_change = 0.1573
    # Carbon macroscopic cross-section = 2820 cm^-1
    macro_scatter_cross_section = Quantity(2820 / units.centimeter)
    Args = namedtuple("Args", ["Na", "Ieff", "Let", "Ss"])
    return Args(Na=atomic_number_density_abs,
        Ieff=resonance_integral,
        Let=average_lethargy_change,
        Ss=macro_scatter_cross_section)


def test_basic_resonance_escape_factor(test_args):
    result = resonance_escape.calculate_resonance_escape_probability(test_args.Na, test_args.Ieff,
        test_args.Let, test_args.Ss)
    assert_approx(result, 0.955)


def test_bad_atomic_number_density(test_args):
    Nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resonance_escape.calculate_resonance_escape_probability(Nb, test_args.Ieff, test_args.Let,
            test_args.Ss)
    with raises(TypeError):
        resonance_escape.calculate_resonance_escape_probability(100, test_args.Ieff, test_args.Let,
            test_args.Ss)
