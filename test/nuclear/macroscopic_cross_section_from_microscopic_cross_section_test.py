from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear import macroscopic_cross_section_from_microscopic_cross_section as macro_cs

# From https://www.nuclear-power.com/nuclear-power/reactor-physics/nuclear-engineering-fundamentals/neutron-nuclear-reactions/macroscopic-cross-section/

Args = namedtuple("Args", ["b", "N"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # carbon microscopic cross-section is 5.01 barn
    microscopic_cross_section = Quantity(5.01e-24 * units.centimeter**2)
    # NC = 2.75*10^22 atoms of carbon/cm^3
    atomic_number_density = Quantity(2.75e22 / units.centimeter**3)
    return Args(b=microscopic_cross_section, N=atomic_number_density)


def test_basic_cross_section(test_args: Args) -> None:
    result = macro_cs.calculate_cross_section(test_args.b, test_args.N)
    assert_equal(result, 0.1377 / units.centimeter)


def test_bad_microscopic_cross_section(test_args: Args) -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        macro_cs.calculate_cross_section(bb, test_args.N)
    with raises(TypeError):
        macro_cs.calculate_cross_section(100, test_args.N)


def test_bad_atomic_number_density(test_args: Args) -> None:
    Nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        macro_cs.calculate_cross_section(test_args.b, Nb)
    with raises(TypeError):
        macro_cs.calculate_cross_section(test_args.b, 100)
