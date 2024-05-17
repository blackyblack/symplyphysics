from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import entropy_derivative_via_volume_derivative as maxwell_relation

# Description
## The isothermal pressure derivative of entropy in a process during which the temperature increased by
## 1 K and the volume dropped by 0.1 m**3 is (dS/dp)_T = 0.1 J/(K*Pa).

Args = namedtuple("Args", "dv dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dv = Quantity(-0.1 * units.meter**3)
    dt = Quantity(1 * units.kelvin)
    return Args(dv=dv, dt=dt)


def test_law(test_args: Args) -> None:
    result = maxwell_relation.calculate_entropy_differential(test_args.dv, test_args.dt)
    assert_equal(result, 0.1 * units.joule / (units.kelvin * units.pascal))


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        maxwell_relation.calculate_entropy_differential(vb, test_args.dt)
    with raises(TypeError):
        maxwell_relation.calculate_entropy_differential(100, test_args.dt)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        maxwell_relation.calculate_entropy_differential(test_args.dv, tb)
    with raises(TypeError):
        maxwell_relation.calculate_entropy_differential(test_args.dv, 100)
