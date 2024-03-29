from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.optics import irradiance_of_light_after_polarizer as malus_law

# Polarized light of unit intensity (I0 = 1 W/m^2) passes through a non-perfect polarizer (k = 0.5) under a 30-degree angle
# in which case the result irradiance is 0.375 W/m^2

Args = namedtuple("Args", ["I0", "k", "phi"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    I0 = Quantity(1 * units.watt / units.meter**2)
    k = 0.5
    phi = Quantity(30 * units.degree)
    return Args(I0=I0, k=k, phi=phi)


def test_basic_irradiance(test_args: Args) -> None:
    result = malus_law.calculate_irradiance(test_args.I0, test_args.k, test_args.phi)
    assert_equal(result, 0.375 * units.watt / units.meter**2)


def test_bad_irradiance_initial(test_args: Args) -> None:
    I0 = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        malus_law.calculate_irradiance(I0, test_args.k, test_args.phi)
    with raises(TypeError):
        malus_law.calculate_irradiance(100, test_args.k, test_args.phi)


def test_bad_transparency_coefficient(test_args: Args) -> None:
    k = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        malus_law.calculate_irradiance(test_args.I0, k, test_args.phi)


def test_bad_polarization_angle(test_args: Args) -> None:
    phi = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        malus_law.calculate_irradiance(test_args.I0, test_args.k, phi)
