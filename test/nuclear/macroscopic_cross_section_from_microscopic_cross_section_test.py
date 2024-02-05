from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.nuclear import macroscopic_cross_section_from_microscopic_cross_section as macro_cs

# From https://www.nuclear-power.com/nuclear-power/reactor-physics/nuclear-engineering-fundamentals/neutron-nuclear-reactions/macroscopic-cross-section/


@fixture(name="test_args")
def test_args_fixture():
    # carbon microscopic cross-section is 5.01 barn
    microscopic_cross_section = Quantity(5.01e-24 * units.centimeter**2)
    # NC = 2.75*10^22 atoms of carbon/cm^3
    atomic_number_density = Quantity(2.75e22 / units.centimeter**3)
    Args = namedtuple("Args", ["b", "N"])
    return Args(b=microscopic_cross_section, N=atomic_number_density)


def test_basic_cross_section(test_args):
    result = macro_cs.calculate_cross_section(test_args.b, test_args.N)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.length)
    result_cross_section = convert_to(result, 1 / units.centimeter).evalf(4)
    assert_approx(result_cross_section, 0.1377)


def test_bad_microscopic_cross_section(test_args):
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        macro_cs.calculate_cross_section(bb, test_args.N)
    with raises(TypeError):
        macro_cs.calculate_cross_section(100, test_args.N)


def test_bad_atomic_number_density(test_args):
    Nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        macro_cs.calculate_cross_section(test_args.b, Nb)
    with raises(TypeError):
        macro_cs.calculate_cross_section(test_args.b, 100)
