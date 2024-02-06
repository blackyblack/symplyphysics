from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_approx,
)
from symplyphysics.laws.chemistry import avogadro_number_from_mole_count


@fixture(name="test_args")
def test_args_fixture():
    mole_count = Quantity(5 * units.mole)
    Args = namedtuple("Args", ["M"])
    return Args(M=mole_count)


def test_basic_particles_count(test_args):
    result = avogadro_number_from_mole_count.calculate_particles_count(test_args.M)
    assert isinstance(result, int)
    assert_approx(result, 3.011e24)


def test_bad_mole_count():
    Mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        avogadro_number_from_mole_count.calculate_particles_count(Mb)
    with raises(TypeError):
        avogadro_number_from_mole_count.calculate_particles_count(100)
