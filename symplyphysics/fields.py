from sympy.vector.operators import _get_coord_systems

def apply_field(field_, trajectory_):
    coord_sys = _get_coord_systems(trajectory_)
    coord_sys = next(iter(coord_sys))
    base_scalars = coord_sys.base_scalars()
    trajectory_projected = trajectory_.to_matrix(coord_sys)
    field_app = field_
    for i, p in enumerate(trajectory_projected):
        field_app = field_app.subs(base_scalars[i], p)
    return field_app