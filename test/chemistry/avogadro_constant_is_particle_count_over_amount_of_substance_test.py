from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.chemistry import avogadro_constant_is_particle_count_over_amount_of_substance

Args = namedtuple("Args", ["M"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mole_count = Quantity(5 * units.mole)
    return Args(M=mole_count)


def test_basic_particles_count(test_args: Args) -> None:
    result = avogadro_constant_is_particle_count_over_amount_of_substance.calculate_particles_count(test_args.M)
    assert isinstance(result, int)
    assert_equal(result, 3.011e24)


def test_bad_mole_count() -> None:
    Mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        avogadro_constant_is_particle_count_over_amount_of_substance.calculate_particles_count(Mb)
    with raises(TypeError):
        avogadro_constant_is_particle_count_over_amount_of_substance.calculate_particles_count(100)
