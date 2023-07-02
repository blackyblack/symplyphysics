from typing import Any, Optional
from sympy.vector import Vector as SymVector, express
from sympy.vector.operators import _get_coord_systems

from ..coordinate_systems.coordinate_systems import CoordinateSystem


# Contains list of SymPy expressions or any numbers as components.
# Contains coordinate system to prevent using vector arithmetics with non compatible
# coordinate systems.
# Vector assumes to have [0, 0, 0] point of origin. All vector related laws work with this assumption
# add do not differentiate between position shifted vectors. But it does not allow to properly rebase
# vectors to another coordinate system.
# Therefore for all physical applications vectors should assume various origin point and should be
# defined dynamically, eg [C.x, C.y] or [parameter1, parameter2].
class Vector:
    _components: list[Any]
    #NOTE: 4 and higher dimensional vectors are not supported cause of using CoordSys3D
    #      to allow rebasing vector coordinate system.
    _coordinate_system: Optional[CoordinateSystem] = None

    def __init__(self, components: list[Any], coordinate_system: Optional[CoordinateSystem] = None):
        self._components = components
        self._coordinate_system = coordinate_system

    @property
    def coordinate_system(self) -> Optional[CoordinateSystem]:
        return self._coordinate_system

    @property
    def components(self) -> list[Any]:
        return self._components


# Converts SymPy Vector to Vector
# SymPy vector is an expression that looks like C.i + C.j, where C is CoordSys3D
def vector_from_sympy_vector(sympy_vector_: SymVector,
    coordinate_system: Optional[CoordinateSystem] = None) -> Vector:
    if sympy_vector_ == SymVector.zero:
        return Vector([], coordinate_system)
    coord_system_set = _get_coord_systems(sympy_vector_)
    coord_system = None
    if len(coord_system_set) > 1:
        coord_sys_names = [str(c) for c in coord_system_set]
        raise TypeError(f"Different coordinate systems in expression: {str(coord_sys_names)}")
    if len(coord_system_set) > 0:
        coord_system = next(iter(coord_system_set))
        if coordinate_system is not None and coord_system != coordinate_system.coord_system:
            raise TypeError(
                f"Different coordinate systems in expression and argument: {str(coord_system)} vs {str(coordinate_system.coord_system)}"
            )
    as_matrix = sympy_vector_.to_matrix(coord_system)
    components = [e for e in as_matrix]
    return Vector(components, coordinate_system)


# Converts Vector to SymPy Vector
def sympy_vector_from_vector(vector_: Vector) -> SymVector:
    result_vector = SymVector.zero
    if vector_.coordinate_system is None:
        return result_vector
    if vector_.coordinate_system.coord_system is None:
        return result_vector
    base_vectors = vector_.coordinate_system.coord_system.base_vectors()
    for idx in range(min(len(base_vectors), len(vector_.components))):
        result_vector = result_vector + base_vectors[idx] * vector_.components[idx]
    return result_vector


# Convert vector coordinate system to new basis and construct new vector.
# Rebased vector should be the same as old vector but in new coordinate system.
def vector_rebase(vector_: Vector, coordinate_system: Optional[CoordinateSystem] = None) -> Vector:
    # Simply set new coordinate system if vector cannot be rebased
    if coordinate_system is None or vector_.coordinate_system is None:
        return Vector(vector_.components, coordinate_system)
    if coordinate_system.coord_system is None or vector_.coordinate_system.coord_system is None:
        return Vector(vector_.components, coordinate_system)
    return _extended_express(vector_, coordinate_system)


def _extended_express(vector_: Vector, system_to: CoordinateSystem) -> Vector:
    # Simply change coordinate system, if nothing to transform from
    if vector_.coordinate_system is None:
        return Vector(vector_.components, system_to)
    if vector_.coordinate_system.coord_system_type != system_to.coord_system_type:
        new_scalars = list(
            vector_.coordinate_system.transformation_to_system(system_to.coord_system_type))
        # now take each component of vector and assign them to base_scalars, eg x, y, z
        # replace each component of new_scalars with assigned x, y, z, eg r -> x, theta -> y
        # build new vector from these components
        for i, scalar in enumerate(vector_.coordinate_system.coord_system.base_scalars()):
            new_component = 0 if i >= len(vector_.components) else vector_.components[i]
            for j, old_scalar in enumerate(new_scalars):
                new_scalars[j] = old_scalar.subs(scalar, new_component)
        vector_ = Vector(new_scalars, vector_.coordinate_system)

    # We do not want to maintain own vector transformation functions, so
    # we convert our vector to SymPy format, transform it and convert back to Vector.
    sympy_vector = sympy_vector_from_vector(vector_)
    transformed_vector_sympy = express(sympy_vector, system_to.coord_system, None, variables=True)
    return vector_from_sympy_vector(transformed_vector_sympy, system_to)
