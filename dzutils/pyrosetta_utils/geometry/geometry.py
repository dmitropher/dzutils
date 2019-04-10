import numpy as _np


def projection(point, plane_point, plane_normal):
    """
    Returns the projection of the given point onto the plane

    vector math cantrips:
    unit_normal dot vector to point is hyp*cos(theta) of a right triangle
    b/w the plane point and the point.

    The orthogonal projection is just point - (normal * height)
    """
    unit_normal = plane_normal / _np.linalg.norm(plane_normal)
    dist = _np.inner(unit_normal, (point - plane_point))
    projection = point - unit_normal * dist
    return projection
