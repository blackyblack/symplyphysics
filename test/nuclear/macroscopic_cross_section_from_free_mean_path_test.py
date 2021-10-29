from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.nuclear import macroscopic_cross_section_from_free_mean_path as macro_cs

@fixture
def test_args():
    mean_free_path = units.Quantity('mean_free_path')
    SI.set_quantity_dimension(mean_free_path, units.length)
    # boron carbide mean free path
    SI.set_quantity_scale_factor(mean_free_path, 0.012 * units.centimeter)

    Args = namedtuple('Args', ['y'])
    return Args(y = mean_free_path)

def test_basic_cross_section(test_args):
    result = macro_cs.calculate_cross_section(test_args.y)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**-1)

    result_cross_section = convert_to(result, units.centimeter**-1).subs(units.centimeter, 1).evalf(2)
    # boron carbide macroscopic cross-section is 84.3 cm^-1
    assert result_cross_section == approx(84.3, 0.1)

def test_bad_microscopic_cross_section():
    yb = units.Quantity('yb')
    SI.set_quantity_dimension(yb, units.time)
    SI.set_quantity_scale_factor(yb, 3 * units.second)

    with raises(errors.UnitsError):
        macro_cs.calculate_cross_section(yb)

    with raises(TypeError):
        macro_cs.calculate_cross_section(100)
