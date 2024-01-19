from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, Quantity, SI, convert_to, errors)
from symplyphysics.laws.electricity import energy_of_an_electron_in_a_hydrogen_atom as energy_law

# Description
## The energy value is known for the radius of the first Bohr orbit of the electron. It is equal to 13.6 [eV].
## https://educon.by/index.php/formuly/phystheory


@fixture(name="test_args")
def test_args_fixture():
    radius_of_electron = Quantity(5.3e-11 * units.meter)
    Args = namedtuple("Args", ["radius_of_electron"])
    return Args(radius_of_electron=radius_of_electron)


def test_basic_energy_of_electron(test_args):
    result = energy_law.calculate_energy_of_electron(test_args.radius_of_electron)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_1 = convert_to(result, units.electronvolt).evalf(5)
    assert result_1 == approx(13.6, 0.01)


def test_bad_radius_of_electron(test_args):
    radius_of_electron = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy_of_electron(radius_of_electron)
    with raises(TypeError):
        energy_law.calculate_energy_of_electron(100)
