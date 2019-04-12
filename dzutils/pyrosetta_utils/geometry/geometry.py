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
    dist = _np.dot((point - plane_point), unit_normal)
    projection = point - dist * unit_normal
    return projection


def closest_point_on_circle(point, center, normal, radius):
    """
    Returns the closest point on the circle given to the search point

    Takes a center, normal, and radius to define the circle
    """
    proj = projection(point, center, normal)
    to_proj = proj - center
    unit = to_proj / _np.linalg.norm(to_proj)
    closest_point = unit * radius + center
    return closest_point
