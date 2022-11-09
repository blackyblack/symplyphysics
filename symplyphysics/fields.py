from typing import List
from sympy.vector import CoordSys3D


# Applies field to a trajectory (surface) - calls field functions with each element of the trajectory as parameter.
# coord_system_ - CoordSys3D coordinate system to work with. Supports Cartesian, Spherical, Cylindrical coordinates.
# field_ - list with functions that map base coordinates to a vector. Each function should have the same number of 
#          parameters as base coordinates in coord_system_. Can contain 0 instead of empty function.
# trajectory_ - list of expressions that correspond to a function in some space. Each element of list
#               corresponds to base coordinate in coord_system_.
# return - list that represents vector parametrized by trajectory parameters.
def apply_field(coord_system_: CoordSys3D, field_: List, trajectory_: List) -> List:
    base_scalars = coord_system_.base_scalars()
    trajectory_size = len(trajectory_)
    field_args = []
    for i in range(len(base_scalars)):
        if i >= trajectory_size:
            field_args.append(0)
            continue
        field_args.append(trajectory_[i])

    field_app = field_.copy()
    for j, f in enumerate(field_app):
        if f == 0: continue
        field_app[j] = f(*field_args)
    return field_app
