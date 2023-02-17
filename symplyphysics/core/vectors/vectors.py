from typing import Any, List
from sympy import Expr, nan
from sympy.vector import CoordSys3D, Vector as SympyVector, express
from sympy.vector.operators import _get_coord_systems


# Detects coordinate system being used in SymPy Vector.
# Throws an exception if multiple coordinate systems are detected.
def extract_coord_system_from_sympy_vector(sympy_vector_: SympyVector) -> CoordSys3D:
    if not isinstance(sympy_vector_, Expr): return None
    coord_system_set = _get_coord_systems(sympy_vector_)
    if len(coord_system_set) > 1:
        coord_sys_names = [str(c) for c in coord_system_set]
        raise TypeError(f"Different coordinate systems in expression: {str(coord_sys_names)}")
    return None if len(coord_system_set) == 0 else next(iter(coord_system_set))

# Contains list of SymPy expressions or any numbers as components.
# Contains coordinate system to prevent using vector arithmetics with non compatible
# coordinate systems.
# Vector assumes to have [0, 0, 0] point of origin. All vector related laws work with this assumption
# add do not differentiate between position shifted vectors. But it does not allow to properly rebase
# vectors to another coordinate system.
# Therefore for all physical applications vectors should assume various origin point and should be
# defined dynamically, eg [C.x, C.y] or [parameter1, parameter2].
class Vector:
    _components: List[Any] = []
    #NOTE: 4 and higher dimensional vectors are not supported cause of using CoordSys3D
    #      to allow vector coordinate system rebasing.
    _coord_system: CoordSys3D = None

    def __init__(self, components: List[Any]=[], coord_system: CoordSys3D=None):
        self._components = components
        self._coord_system = coord_system

    @property
    def coord_system(self) -> CoordSys3D:
        return self._coord_system

    @property
    def components(self):
        return self._components

# Converts SymPy Vector to Vector
# SymPy vector is an expression that looks like C.i + C.j, where C is CoordSys3D
def vector_from_sympy_vector(sympy_vector_: Any) -> Vector:
    if sympy_vector_ == SympyVector.zero: return Vector([])
    coord_system = extract_coord_system_from_sympy_vector(sympy_vector_)
    # Do not convert if cannot detect SymPy Vector class
    if coord_system is None or not sympy_vector_.is_Vector: return Vector([sympy_vector_])
    as_matrix = sympy_vector_.to_matrix(coord_system)
    components = []
    for e in as_matrix:
        components.append(e)
    return Vector(components, coord_system)

# Converts Vector to SymPy Vector
def sympy_vector_from_vector(vector_: Vector) -> SympyVector:
    result_vector = SympyVector.zero
    if vector_.coord_system is None: return result_vector
    base_vectors = vector_.coord_system.base_vectors()
    for idx in range(min(len(base_vectors), len(vector_.components))):
        result_vector = result_vector + base_vectors[idx] * vector_.components[idx]
    return result_vector

# Convert vector coordinate system to new basis and construct new vector.
# Rebased vector should be the same as old vector but in new coordinate system.
def vector_rebase(vector_: Vector, coord_system: CoordSys3D=None) -> Vector:
    # Simply set new coordinate system if vector cannot be rebased
    if coord_system is None or vector_.coord_system is None:
        return Vector(vector_.components, coord_system)

    # We do not want to maintain own vector transformation functions, so
    # we convert our vector to SymPy format, transform it and convert back to Vector.
    sympy_vector = sympy_vector_from_vector(vector_)
    transformed_vector = express(sympy_vector, coord_system, variables=True)
    return vector_from_sympy_vector(transformed_vector)

def extended_express(expr: Expr, system: CoordSys3D, system2: CoordSys3D=None, variables=False):
    if system._transformation_from_parent_lambda is not None:
        return _extended_express(expr, system)
    return express(expr, system, system2, variables)

#TODO: this transformation does not properly work for coordinate system translations and rotations. Should be
#      rewritten to support both standard express functionality and non-cartesian coordinate system
#      transformations.
# This is experimental feature from SymPy (not merged yet)
def _extended_express(expr: Expr, coordsys: CoordSys3D):
    """
    express for vectors allowing transformation between
    cartesian (x,y,z) and non-cartesian frames (eg cylindrical)
    Parameters
    ==========
    expr : Vector
        The expression to re-express in CoordSys3D 'coordsys'
    system: CoordSys3D
        The coordinate system the expr is to be expressed in
    Examples
    =========
    >>> from sympy.vector import extended_express, CoordSys3D
    >>> from sympy import pi
    >>> N = CoordSys3D('N')
    >>> C = N.create_new('C', transformation='cylindrical', \
        vector_names = ('r', 'theta', 'z'))
    >>> extended_express(5*N.i,C)
    5*C.r
    >>> extended_express(5*N.j,C)
    5*C.r + pi/2*C.theta
    >>> extended_express(10*N.k,C)
    10*C.z
    >>> extended_express(6*C.r,N)
    6*N.i
    >>> extended_express(6*C.r + pi/2*C.theta,N)
    6*N.j
    Where a coordinate is undefined (eg theta with
    extended_express(N.k), it is mapped to zero.
    The non-cartesian frame must be defined as a child of the
    cartesian frame as above.
    The output is expressed in terms of the base vectors of the
    derived coordinate systems. Currently they default to i,j
    and k which are confusing for non-cartesian systems so it
    is recommended to explicitly name them as above.
    See Also
    ========
    CoordSys3D.transformation_from_parent,
    CoordSys3D.transformation_to_parent
    """
    if not isinstance(coordsys, CoordSys3D):
        raise TypeError("system should be a CoordSys3D instance")

    parts = expr.separate()

    # get different coordinate frames in vector, cry if more than one
    if isinstance(parts, dict):
        # we have been given zero vector
        if(len(parts) == 0): return SympyVector.zero
        if(len(parts) > 1):  raise ValueError("Does not support \
            expressions containing multiple base coordinate frames")

        foundframe = tuple(parts.keys())[0]
        foundvector = tuple(parts.values())[0]
    else:
        foundframe = parts.system
        foundvector = parts

    # we should now have found the only frame in the vector,
    # and we now have to convert it to the given frame

    from_frame = foundframe
    to_frame = coordsys

    from_coeffs1 = from_frame.base_vectors()
    from_coeffs2 = from_frame.base_scalars()
    to_coeffs = to_frame.base_vectors()    # output with vectors

    if from_frame._parent == to_frame:
        transform_function = from_frame._transformation_lambda
    else:
        if to_frame._parent == from_frame:
            transform_function =  \
                to_frame._transformation_from_parent_lambda
        else:
            raise ValueError("Cannot link Coordinate frames")
    args = []

    for i, j in zip(from_coeffs1, from_coeffs2):
        # could be expressed in either for inertial frame
        coeff1 = foundvector.coeff(i)
        # understand both N.i and N.x (or C.r and C.i)
        coeff2 = foundvector.coeff(j)
        if coeff1 != 0 and coeff2 != 0:
            raise ValueError("Cannot express vector with both base \
                    vectors and base scalars - check your vector")
        args.append(coeff1 + coeff2) # only one can be non-zero

    vals = transform_function(*args)

    ans = SympyVector.zero
    for v, c in zip(vals, to_coeffs):
        if v is nan:  # v is nan ie infinity from an arctan
            v = 0
        # use to solve for cylindrical coords
        # where theta is undefined
        ans += v * c
    return ans