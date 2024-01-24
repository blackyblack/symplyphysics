from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
    prefixes,
)
from symplyphysics.laws.thermodynamics import energy_to_vaporization_from_mass as amount_energy

# https://easyfizika.ru/zadachi/termodinamika/na-zazhzhennuyu-spirtovku-s-kpd-60-postavili-sosud-s-500-g-vody-pri-20-c-cherez-kakoe/
# If mass of water equal 20 g and specific heat of vaporization equal 29 MJ/kg,
# the energy to vaporization should be equal to 0.58 MJ


@fixture(name="test_args")
def test_args_fixture():
    k_v = Quantity(29 * prefixes.mega * units.joule / units.kilogram)
    m = Quantity(20 * units.grams)
    Args = namedtuple("Args", ["k_v", "m"])
    return Args(k_v=k_v, m=m)


def test_basic_amount(test_args):
    result = amount_energy.calculate_amount_energy(test_args.k_v, test_args.m)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, prefixes.mega * units.joules).evalf(5)
    assert result_energy == approx(0.58, 0.01)


def test_bad_specific_heat_vaporization(test_args):
    k_vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(k_vb, test_args.m)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(100, test_args.m)


def test_bad_body_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.k_v, mb)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.k_v, 100)
