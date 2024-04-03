from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import thermodynamic_compressibility as compressibility_def

# Description
## During a process the volume of a body changed from 1 m**3 to 0.999 m**3 while the pressure
## in the body changed from 1 Pa to 1.005 Pa. The value of the thermodynamic compressibility
## coefficient at the latter pressure-volume point amounts to about 0.2 1/Pa.

Args = namedtuple("Args", "v1 v2 p1 p2")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v1 = Quantity(1 * units.meter**3)
    v2 = Quantity(0.999 * units.meter**3)
    p1 = Quantity(1 * units.pascal)
    p2 = Quantity(1.005 * units.pascal)
    return Args(v1=v1, v2=v2, p1=p1, p2=p2)


def test_definition(test_args: Args) -> None:
    result = compressibility_def.calculate_compressibility(test_args.v1, test_args.v2, test_args.p1, test_args.p2)
    assert_equal(result, 0.2 / units.pascal, tolerance=2e-3)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        compressibility_def.calculate_compressibility(vb, test_args.v2, test_args.p1, test_args.p2)
    with raises(TypeError):
        compressibility_def.calculate_compressibility(100, test_args.v2, test_args.p1, test_args.p2)
    with raises(errors.UnitsError):
        compressibility_def.calculate_compressibility(test_args.v1, vb, test_args.p1, test_args.p2)
    with raises(TypeError):
        compressibility_def.calculate_compressibility(test_args.v1, 100, test_args.p1, test_args.p2)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        compressibility_def.calculate_compressibility(test_args.v1, test_args.v2, pb, test_args.p2)
    with raises(TypeError):
        compressibility_def.calculate_compressibility(test_args.v1, test_args.v2, 100, test_args.p2)
    with raises(errors.UnitsError):
        compressibility_def.calculate_compressibility(test_args.v1, test_args.v2, test_args.p1, pb)
    with raises(TypeError):
        compressibility_def.calculate_compressibility(test_args.v1, test_args.v2, test_args.p1, 100)
