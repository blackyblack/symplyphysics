from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import (
    elastic_energy_density_of_compression_via_strain as law,)

Args = namedtuple("Args", "e s")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = Quantity(10 * units.pascal)
    s = 0.05
    return Args(e=e, s=s)


def test_law(test_args: Args) -> None:
    result = law.calculate_elastic_energy_density(test_args.e, test_args.s)
    assert_equal(result, 13e-3 * units.joule / units.meter**3, relative_tolerance=4e-2)


def test_bad_pressure(test_args: Args) -> None:
    eb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_elastic_energy_density(eb, test_args.s)
    with raises(TypeError):
        law.calculate_elastic_energy_density(100, test_args.s)


def test_bad_strain(test_args: Args) -> None:
    eb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_elastic_energy_density(test_args.e, eb)
