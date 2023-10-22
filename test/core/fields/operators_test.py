from collections import namedtuple
from pytest import fixture
from sympy import cos, exp, sin, Symbol as SymSymbol
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.fields.operators import divergence_operator


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    parameter1 = SymSymbol("parameter1")
    parameter2 = SymSymbol("parameter2")
    Args = namedtuple("Args", ["C", "parameter1", "parameter2"])
    return Args(C=C, parameter1=parameter1, parameter2=parameter2)


def test_basic_divergence(test_args):
    field = VectorField(
        lambda point: [exp(point.x) * cos(point.y),
        exp(point.x) * sin(point.y), point.z], test_args.C)
    result = divergence_operator(field)
    x = field.coordinate_system.coord_system.base_scalars()[0]
    y = field.coordinate_system.coord_system.base_scalars()[1]
    assert result == 2 * exp(x) * cos(y) + 1


#TODO: divergence_operator tests
#TODO: curl_operator tests
