from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.waveguides import attenuation_coefficient_in_metal as coefficient_law

# Description
## Parameters of the coaxial waveguide: the diameter of the inner wire is 1.35 millimeters, the diameter of the outer conductor
## is 9.0 millimeters, the relative permeability of the dielectric is 1, the relative permittivity of the dielectric is 2.2.
## The surface resistance of the outer conductor and the surface resistance of the inner conductor are 2.576 milliohm.
## The attenuation coefficient will be 0.0013 [1 / meter].
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", [
    "relative_permittivity", "relative_permeability", "surface_resistance_outer", "surface_resistance_inner", "outer_diameter",
    "inner_diameter"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 2.2
    relative_permeability = 1
    surface_resistance_outer = Quantity(2.576 * prefixes.milli * units.ohm)
    surface_resistance_inner = Quantity(2.576 * prefixes.milli * units.ohm)
    outer_diameter = Quantity(9 * units.millimeter)
    inner_diameter = Quantity(1.35 * units.millimeter)

    return Args(relative_permittivity=relative_permittivity,
        relative_permeability=relative_permeability,
        surface_resistance_outer=surface_resistance_outer,
        surface_resistance_inner=surface_resistance_inner,
        outer_diameter=outer_diameter,
        inner_diameter=inner_diameter)


def test_basic_attenuation_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity,
        test_args.relative_permeability,
        test_args.surface_resistance_outer, test_args.surface_resistance_inner, test_args.outer_diameter,
        test_args.inner_diameter)
    assert_equal(result, 0.0013 * (1 / units.meter))


def test_bad_relative_permittivity(test_args: Args) -> None:
    relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(relative_permittivity, test_args.relative_permeability,
            test_args.surface_resistance_outer, test_args.surface_resistance_inner, test_args.outer_diameter,
            test_args.inner_diameter)


def test_bad_relative_permeability(test_args: Args) -> None:
    relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity, relative_permeability,
            test_args.surface_resistance_outer, test_args.surface_resistance_inner, test_args.outer_diameter,
            test_args.inner_diameter)


def test_bad_surface_resistance(test_args: Args) -> None:
    bad_surface_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity, test_args.relative_permeability,
            bad_surface_resistance, test_args.surface_resistance_inner, test_args.outer_diameter,
            test_args.inner_diameter)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity, test_args.relative_permeability, 100,
            test_args.surface_resistance_inner, test_args.outer_diameter, test_args.inner_diameter)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity, test_args.relative_permeability,
            test_args.surface_resistance_outer, bad_surface_resistance, test_args.outer_diameter,
            test_args.inner_diameter)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity, test_args.relative_permeability,
            test_args.surface_resistance_outer, 100, test_args.outer_diameter,
            test_args.inner_diameter)


def test_bad_radius(test_args: Args) -> None:
    bad_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity, test_args.relative_permeability, test_args.surface_resistance_outer, test_args.surface_resistance_inner, bad_radius, test_args.inner_diameter)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity, test_args.relative_permeability, test_args.surface_resistance_outer, test_args.surface_resistance_inner, 100, test_args.inner_diameter)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity, test_args.relative_permeability, test_args.surface_resistance_outer, test_args.surface_resistance_inner, test_args.outer_diameter, bad_radius)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity, test_args.relative_permeability, test_args.surface_resistance_outer, test_args.surface_resistance_inner, test_args.outer_diameter, 100)
    with raises(ValueError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity, test_args.relative_permeability, test_args.surface_resistance_outer, test_args.surface_resistance_inner, test_args.inner_diameter, test_args.outer_diameter)
