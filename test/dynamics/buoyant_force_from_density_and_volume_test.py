from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.dynamics import buoyant_force_from_density_and_volume as archimedes_law

Args = namedtuple("Args", ["V", "pf"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # water density
    pf = Quantity(1000 * units.kilogram / units.meter**3)
    V = Quantity(0.2 * units.meter**3)
    return Args(V=V, pf=pf)


def test_basic_force(test_args: Args) -> None:
    result = archimedes_law.calculate_force_buoyant(test_args.pf, test_args.V)
    assert_equal(result, 1961.3 * units.newton)


def test_bad_density(test_args: Args) -> None:
    pb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        archimedes_law.calculate_force_buoyant(pb, test_args.V)
    with raises(TypeError):
        archimedes_law.calculate_force_buoyant(100, test_args.V)


def test_bad_volume(test_args: Args) -> None:
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        archimedes_law.calculate_force_buoyant(test_args.pf, Vb)
    with raises(TypeError):
        archimedes_law.calculate_force_buoyant(test_args.pf, 100)
