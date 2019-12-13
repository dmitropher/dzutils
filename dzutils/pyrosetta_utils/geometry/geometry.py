import numpy as _np


def rmsd(atom_coords_1, atom_coords_2):
    """
    Takes two coord arrays of dim 2 or more and gets the rmsd

    All numpy arrays, broadcastable, etc
    """
    return _np.sqrt(
        _np.sum(
            _np.square(
                _np.linalg.norm(
                    atom_coords_2 - atom_coords_1,
                    axis=len(atom_coords_2.shape) - 2,
                )
            )
            # shape -2 index here should be the number of points
            / atom_coords_2.shape[-2],
            axis=len(atom_coords_2.shape) - 2,
        )
    )


def unit_vec(vec):
    """
    Returns the vector divided by its mag
    """
    return vec / _np.linalg.norm(vec)


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


def vector_theta(vec_1, vec_2, radians=False):
    """
    Returns the theta between vec_1 and vec_2
    """
    cos_theta = _np.dot(unit_vec(vec_1), unit_vec(vec_2))
    if radians:
        return _np.arccos(_np.clip(cos_theta, -1, 1))

    return _np.degrees(_np.arccos(_np.clip(cos_theta, -1, 1)))
