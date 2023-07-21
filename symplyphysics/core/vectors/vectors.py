from typing import Sequence, TypeAlias
from sympy import Expr, Add, Mul
from sympy.vector import VectorAdd, VectorMul, Vector as SymVector, express
from sympy.vector.operators import _get_coord_systems

from ...core.symbols.quantities import Quantity
from ..coordinate_systems.coordinate_systems import CoordinateSystem

VectorComponent: TypeAlias = Expr | float | Quantity


# Contains list of SymPy expressions or any numbers as components.
# Contains coordinate system to prevent using vector arithmetics with non compatible
# coordinate systems.
# Vector assumes to have [0, 0, 0] point of origin. All vector related laws work with this assumption
# add do not differentiate between position shifted vectors. But it does not allow to properly rebase
# vectors to another coordinate system.
# Therefore for all physical applications vectors should assume various origin point and should be
# defined dynamically, eg [C.x, C.y] or [parameter1, parameter2].
#NOTE: use Python lists when coordinate system does not matter
class Vector:
    #NOTE: 4 and higher dimensional vectors are not supported cause of using CoordSys3D
    #      to allow rebasing vector coordinate system.
    _coordinate_system: CoordinateSystem
    _components: list[VectorComponent]

    def __init__(self, coordinate_system: CoordinateSystem, components: Sequence[VectorComponent]):
        self._coordinate_system = coordinate_system
        self._components = list(components)

    @property
    def coordinate_system(self) -> CoordinateSystem:
        return self._coordinate_system

    @property
    def components(self) -> Sequence[VectorComponent]:
        return self._components


# Converts SymPy Vector to Vector
# SymPy vector is an expression that looks like C.i + C.j, where C is CoordSys3D
def vector_from_sympy_vector(sympy_vector_: SymVector,
    coordinate_system: CoordinateSystem) -> Vector:
    if sympy_vector_ == SymVector.zero:
        return Vector(coordinate_system, [])
    coord_system_set = _get_coord_systems(sympy_vector_)
    coord_system = None
    if len(coord_system_set) > 1:
        coord_sys_names = [str(c) for c in coord_system_set]
        raise TypeError(f"Different coordinate systems in expression: {str(coord_sys_names)}")
    if len(coord_system_set) > 0:
        coord_system = next(iter(coord_system_set))
        if coord_system != coordinate_system.coord_system:
            raise TypeError(
                f"Different coordinate systems in expression and argument: {str(coord_system)} vs {str(coordinate_system.coord_system)}"
            )
    components = list(sympy_vector_.to_matrix(coord_system))
    return Vector(coordinate_system, components)


# Converts Vector to SymPy Vector
def sympy_vector_from_vector(vector_: Vector) -> SymVector:
    result_vector = SymVector.zero
    base_vectors = vector_.coordinate_system.coord_system.base_vectors()
    for idx in range(min(len(base_vectors), len(vector_.components))):
        result_vector = result_vector + base_vectors[idx] * vector_.components[idx]
    return result_vector


# Convert vector coordinate system to new basis and construct new vector.
# Rebased vector should be the same as old vector but in new coordinate system.
def vector_rebase(vector_: Vector, coordinate_system: CoordinateSystem) -> Vector:
    if vector_.coordinate_system.coord_system_type != coordinate_system.coord_system_type:
        new_scalars = list(
            vector_.coordinate_system.transformation_to_system(coordinate_system.coord_system_type))
        # now take each component of vector and assign them to base_scalars, eg x, y, z
        # replace each component of new_scalars with assigned x, y, z, eg r -> x, theta -> y
        # build new vector from these components
        for i, scalar in enumerate(vector_.coordinate_system.coord_system.base_scalars()):
            new_component = 0 if i >= len(vector_.components) else vector_.components[i]
            for j, old_scalar in enumerate(new_scalars):
                new_scalars[j] = old_scalar.subs(scalar, new_component)
        vector_ = Vector(vector_.coordinate_system, new_scalars)

    # We do not want to maintain own vector transformation functions, so
    # we convert our vector to SymPy format, transform it and convert back to Vector.
    sympy_vector = sympy_vector_from_vector(vector_)
    transformed_vector_sympy = express(sympy_vector,
        coordinate_system.coord_system,
        None,
        variables=True)
    return vector_from_sympy_vector(transformed_vector_sympy, coordinate_system)


def expr_to_vector(expr: Expr, coordinate_system: CoordinateSystem) -> Vector:
    if isinstance(expr, Mul):
        expr = VectorMul(*expr.args)
    if isinstance(expr, Add):
        expr = VectorAdd(*expr.args)
    if not isinstance(expr, SymVector):
        raise TypeError(f"Expression cannot be converted to SymPy Vector: {str(expr)}")
    vector = vector_from_sympy_vector(expr, coordinate_system)
    components: list[Expr | float | Quantity] = []
    for c in vector.components:
        components.append(Quantity(c))
    return Vector(coordinate_system, components)
