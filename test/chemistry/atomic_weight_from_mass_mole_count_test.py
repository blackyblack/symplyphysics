from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.chemistry import atomic_weight_from_mass_mole_count


@fixture(name="test_args")
def test_args_fixture():
    # molar mass of water is 18.0153 gram / mole
    mass = Quantity(18.0153 * units.gram)
    mole_count = Quantity(1 * units.mole)
    Args = namedtuple("Args", ["m", "N"])
    return Args(m=mass, N=mole_count)


def test_basic_atomic_weight(test_args):
    result = atomic_weight_from_mass_mole_count.calculate_atomic_weight(test_args.m, test_args.N)
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
        units.mass / units.amount_of_substance)
    result_mass = convert_to(result, units.gram / units.mole).evalf(4)
    assert result_mass == approx(18.0153, 0.01)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(mb, test_args.N)
    with raises(TypeError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(100, test_args.N)


def test_bad_mole_count(test_args):
    Nb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(test_args.m, Nb)
    with raises(TypeError):
        atomic_weight_from_mass_mole_count.calculate_atomic_weight(test_args.m, 100)
