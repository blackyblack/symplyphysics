from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.definitions import (
    compressibility_factor_is_deviation_from_ideal_gas as deviation_def,)

# Description
## The compressibility factor of a gas with pressure p = 5.5 kPa, temperature T = 300 K,
## of an amount n = 2 mol, confined in a space of volume V = 1 m**3, is Z = 1.1.

Args = namedtuple("Args", "p v n t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(5.5 * prefixes.kilo * units.pascal)
    v = Quantity(1 * units.meter**3)
    n = Quantity(2 * units.mole)
    t = Quantity(300 * units.kelvin)
    return Args(p=p, v=v, n=n, t=t)


def test_definition(test_args: Args) -> None:
    result = deviation_def.calculate_compressibility_factor(test_args.p, test_args.v, test_args.n,
        test_args.t)
    assert_equal(result, 1.1, relative_tolerance=3e-3)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        deviation_def.calculate_compressibility_factor(pb, test_args.v, test_args.n, test_args.t)
    with raises(TypeError):
        deviation_def.calculate_compressibility_factor(100, test_args.v, test_args.n, test_args.t)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        deviation_def.calculate_compressibility_factor(test_args.p, vb, test_args.n, test_args.t)
    with raises(TypeError):
        deviation_def.calculate_compressibility_factor(test_args.p, 100, test_args.n, test_args.t)


def test_bad_amount_of_substance(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        deviation_def.calculate_compressibility_factor(test_args.p, test_args.v, nb, test_args.t)
    with raises(TypeError):
        deviation_def.calculate_compressibility_factor(test_args.p, test_args.v, 100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        deviation_def.calculate_compressibility_factor(test_args.p, test_args.v, test_args.n, tb)
    with raises(TypeError):
        deviation_def.calculate_compressibility_factor(test_args.p, test_args.v, test_args.n, 100)
