from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.circuits import specific_inductance_of_coaxial_waveguide as inductance_law

## Parameters of the coaxial waveguide: the radius of the inner wire is 1.35 millimeters, the radius of the outer conductor
## is 9.0 millimeters, the relative permeability of the dielectric is 1. The specific inductance will be 379.424 nanohenry per meter.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", ["relative_permeability", "outer_radius", "inner_radius"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permeability = 1
    outer_radius = Quantity(9 * units.millimeter)
    inner_radius = Quantity(1.35 * units.millimeter)
    return Args(relative_permeability=relative_permeability, outer_radius=outer_radius, inner_radius=inner_radius)


def test_basic_specific_inductance(test_args: Args) -> None:
    result = inductance_law.calculate_specific_inductance(test_args.relative_permeability, test_args.outer_radius,
        test_args.inner_radius)
    assert_equal(result, 379.424 * prefixes.nano * units.henry / units.meter)


def test_bad_relative_permeability(test_args: Args) -> None:
    bad_relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inductance_law.calculate_specific_inductance(bad_relative_permeability, test_args.outer_radius, test_args.inner_radius)


def test_bad_radius(test_args: Args) -> None:
    bad_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inductance_law.calculate_specific_inductance(test_args.relative_permeability, bad_radius, test_args.inner_radius)
    with raises(TypeError):
        inductance_law.calculate_specific_inductance(test_args.relative_permeability, 100, test_args.inner_radius)
    with raises(errors.UnitsError):
        inductance_law.calculate_specific_inductance(test_args.relative_permeability, test_args.outer_radius, bad_radius)
    with raises(TypeError):
        inductance_law.calculate_specific_inductance(test_args.relative_permeability, test_args.outer_radius, 100)
