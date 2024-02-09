from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear import diffusion_area_from_diffusion_coefficient as diffusion_area

Args = namedtuple("Args", ["Sa", "D"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # water macroscopic absorption cross-section is 0.022 cm^-1
    macro_abs_cross_section = Quantity(0.022 / units.centimeter)
    # water diffusion coefficient is 0.142 cm
    diffusion_coefficient = Quantity(0.142 * units.centimeter)
    return Args(Sa=macro_abs_cross_section, D=diffusion_coefficient)


def test_basic_diffusion_length(test_args: Args) -> None:
    result = diffusion_area.calculate_diffusion_area(test_args.D, test_args.Sa)
    assert_equal(result, 2.54**2 * units.centimeter**2)


def test_bad_diffusion_coefficient(test_args: Args) -> None:
    Db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_area.calculate_diffusion_area(Db, test_args.Sa)
    with raises(TypeError):
        diffusion_area.calculate_diffusion_area(100, test_args.Sa)


def test_bad_macroscopic_cross_section(test_args: Args) -> None:
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_area.calculate_diffusion_area(test_args.D, Sb)
    with raises(TypeError):
        diffusion_area.calculate_diffusion_area(test_args.D, 100)
