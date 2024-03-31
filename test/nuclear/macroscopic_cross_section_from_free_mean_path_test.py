from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear import macroscopic_cross_section_from_free_mean_path as macro_cs

# From https://www.nuclear-power.com/nuclear-power/reactor-physics/nuclear-engineering-fundamentals/neutron-nuclear-reactions/macroscopic-cross-section/

Args = namedtuple("Args", ["y"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # boron carbide mean free path
    mean_free_path = Quantity(0.01186 * units.centimeter)
    return Args(y=mean_free_path)


def test_basic_cross_section(test_args: Args) -> None:
    result = macro_cs.calculate_cross_section(test_args.y)
    # boron carbide macroscopic cross-section is 84.3 cm^-1
    assert_equal(result, 84.3 / units.centimeter)


def test_bad_microscopic_cross_section() -> None:
    yb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        macro_cs.calculate_cross_section(yb)
    with raises(TypeError):
        macro_cs.calculate_cross_section(100)
