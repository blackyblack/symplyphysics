from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity.circuits import admittance_of_parallel_dipoles as parallel_law

# Description
## Assert we have two resistors with 1/2 Siemens and 1/4 Siemens conductance.
## Accordind to calculator (https://www.chipdip.ru/calc/parallel-resistors) resulting conductance should be 0.75 Siemens.
## Hence resistors are dipoles, test fixture for parallel resistors law should be suitable for parallel dipoles law as well.


@fixture(name="test_args")
def test_args_fixture():
    S1 = Quantity(1 / 2 * units.siemens)
    S2 = Quantity(1 / 4 * units.siemens)
    Args = namedtuple("Args", ["S1", "S2"])
    return Args(S1=S1, S2=S2)


def test_basic_admittance(test_args):
    result = parallel_law.calculate_parallel_admittance([test_args.S1, test_args.S2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.conductance)
    result_admittance = convert_to(result, units.siemens).evalf(3)
    assert result_admittance == approx(0.75, 0.001)


def test_three_resistors_array(test_args):
    S3 = Quantity(1 * units.siemens)
    result = parallel_law.calculate_parallel_admittance([test_args.S1, test_args.S2, S3])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.conductance)
    result_admittance = convert_to(result, units.siemens).evalf(3)
    assert result_admittance == approx(1.75, 0.01)


def test_bad_admittance(test_args):
    Sb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        parallel_law.calculate_parallel_admittance([Sb, test_args.S2])
    with raises(TypeError):
        parallel_law.calculate_parallel_admittance([100, test_args.S2])
    with raises(errors.UnitsError):
        parallel_law.calculate_parallel_admittance([test_args.S1, Sb])
    with raises(TypeError):
        parallel_law.calculate_parallel_admittance([test_args.S1, 100])
    with raises(errors.UnitsError):
        parallel_law.calculate_parallel_admittance([Sb, Sb])
    with raises(TypeError):
        parallel_law.calculate_parallel_admittance([100, 100])
    with raises(TypeError):
        parallel_law.calculate_parallel_admittance(test_args.S1)
