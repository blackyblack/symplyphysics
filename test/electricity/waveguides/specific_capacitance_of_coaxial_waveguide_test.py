from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.waveguides import specific_capacitance_of_coaxial_waveguide as capacitance_law

## Parameters of the coaxial waveguide: the radius of the inner wire is 1.35 millimeters, the radius of the outer conductor
## is 9.0 millimeters, the relative permittivity of the dielectric is 2.2. The specific capacitance will be 64.5 picofarad per meter.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", ["relative_permittivity", "outer_radius", "inner_radius"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 2.2
    outer_radius = Quantity(9 * units.millimeter)
    inner_radius = Quantity(1.35 * units.millimeter)
    return Args(relative_permittivity=relative_permittivity, outer_radius=outer_radius, inner_radius=inner_radius)


def test_basic_specific_capacitance(test_args: Args) -> None:
    result = capacitance_law.calculate_specific_capacitance(test_args.relative_permittivity, test_args.outer_radius,
        test_args.inner_radius)
    assert_equal(result, 64.5 * prefixes.pico * units.farad / units.meter)


def test_bad_relative_permittivity(test_args: Args) -> None:
    bad_relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitance_law.calculate_specific_capacitance(bad_relative_permittivity, test_args.outer_radius, test_args.inner_radius)


def test_bad_radius(test_args: Args) -> None:
    bad_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitance_law.calculate_specific_capacitance(test_args.relative_permittivity, bad_radius, test_args.inner_radius)
    with raises(TypeError):
        capacitance_law.calculate_specific_capacitance(test_args.relative_permittivity, 100, test_args.inner_radius)
    with raises(errors.UnitsError):
        capacitance_law.calculate_specific_capacitance(test_args.relative_permittivity, test_args.outer_radius, bad_radius)
    with raises(TypeError):
        capacitance_law.calculate_specific_capacitance(test_args.relative_permittivity, test_args.outer_radius, 100)
    with raises(ValueError):
        capacitance_law.calculate_specific_capacitance(test_args.relative_permittivity, test_args.inner_radius, test_args.outer_radius)