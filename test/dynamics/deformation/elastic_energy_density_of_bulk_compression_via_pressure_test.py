from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import (
    elastic_energy_density_of_bulk_compression_via_pressure as law,)

Args = namedtuple("Args", "p k")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(5 * prefixes.mega * units.pascal)
    k = Quantity(2 * prefixes.giga * units.pascal)
    return Args(p=p, k=k)


def test_law(test_args: Args) -> None:
    result = law.calculate_elastic_energy_density(test_args.p, test_args.k)
    assert_equal(result, 6.2 * prefixes.kilo * units.joule / units.meter**3, tolerance=8e-3)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_elastic_energy_density(pb, test_args.k)
    with raises(TypeError):
        law.calculate_elastic_energy_density(100, test_args.k)
    with raises(errors.UnitsError):
        law.calculate_elastic_energy_density(test_args.p, pb)
    with raises(TypeError):
        law.calculate_elastic_energy_density(test_args.p, 100)
