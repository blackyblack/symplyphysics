from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.nuclear import macroscopic_cross_section_from_microscopic_cross_section as macro_cs


@fixture
def test_args():
    # carbon microscopic cross-section is 5.01 barn
    microscopic_cross_section = Quantity(5.01e-24 * units.centimeter**2)
    # NC = 1 * 2.75*10^22 atoms of carbon/cm^3
    atomic_number_density = Quantity(2.75e22 / units.centimeter**3)
    Args = namedtuple("Args", ["b", "N"])
    return Args(b=microscopic_cross_section, N=atomic_number_density)


def test_basic_cross_section(test_args):
    result = macro_cs.calculate_cross_section(test_args.b, test_args.N)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**-1)
    result_cross_section = convert_to(result, units.centimeter**-1).subs(units.centimeter,
        1).evalf(2)
    assert result_cross_section == approx(0.14, 0.1)


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
