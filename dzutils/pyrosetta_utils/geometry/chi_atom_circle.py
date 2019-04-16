import numpy as _np
from dzutils.pyrosetta_utils.geometry.homog import (
    np_rot_trans_to_homog as _to_homog,
)
from dzutils.pyrosetta_utils.geometry import (
    closest_point_on_circle as _closest,
    vector_theta as vector_theta,
    unit_vec as unit_vec,
)

# Do for every dock
# save the inverse rotamer circle as just radius, point and unit normal vector
# compute rotation translation to origin (change of basis)
class ChiAtomCircle:
    """
    Describes the circle formed when tracking atom rotation about the given chi

    Works lazily, use getters and setters.
    Will also break lazily if you give it stuff with no chis.
    Decides for you whether you're looking for forward or inverse rotamers
    based on atom numbers (it's not very smart).
    """

    def __init__(self, residue, chi_num, atom="CA"):
        """
        Takes a rosetta residue object, a chi_num and an atom name str or index

        Default atom is CA, best used to search for inverse rotamers
        """
        # set chi,atom,residue variables

        self._chi = chi_num
        self._residue = residue.clone()
        self._chi_atoms = [
            a for a in self._residue.type().chi_atoms()[self._chi]
        ]
        self._atom_xyz = self._residue.xyz(atom)
        self._atom_num = (
            atom
            if isinstance(atom, int)
            else [
                i
                for i, a in enumerate(self._residue.atoms(), 1)
                if a.mm_type() == self._residue.atom(atom).mm_type()
            ][0]
        )
        self._atom_name = (
            atom if isinstance(atom, str) else self._residue.atom_name(atom)
        )
        self._atom = self._residue.atom(atom)

    def compute_circle(self):
        """
        Computes the circle from the residue, chi, and atoms

        Basically draws a right triangle with the vector the target atom as
        the hypotenuse.

        Trig identities to the find the orgin/radius, normal is just the
        unit vector down the chi in either the for or rev rotamer direction.
        """
        # Figure out whether to compute inverse or forward rotamer circle
        self._inverse_rotamer = (
            True if self._atom_num < self._chi_atoms[1] else False
        )
        vec_atoms = (
            (self._chi_atoms[2], self._chi_atoms[1])
            if self._inverse_rotamer
            else (self._chi_atoms[1], self._chi_atoms[2])
        )

        chi_origin = _np.array(
            [coor for coor in self._residue.xyz(vec_atoms[0])]
        )
        atom_xyz = _np.array([coor for coor in self._atom.xyz()])
        chi_vec = (
            _np.array([coor for coor in self._residue.xyz(vec_atoms[1])])
            - chi_origin
        )
        unit_chi_vec = chi_vec / _np.linalg.norm(chi_vec)
        hyp_vec = atom_xyz - chi_origin
        hyp_norm = _np.linalg.norm(hyp_vec)
        unit_hyp_vec = hyp_vec / hyp_norm
        cos_theta = _np.dot(unit_chi_vec, unit_hyp_vec)
        self._normal_vector = unit_chi_vec
        self._circle_origin = chi_origin + unit_chi_vec * cos_theta * hyp_norm
        self._circle_radius = (
            _np.linalg.norm(_np.cross(unit_hyp_vec, self._normal_vector))
            * hyp_norm
        )

    def compute_rt_to_xy(self):
        """
        Returns the homog rt xform to the xy plane

        Useful for projecting another point onto the plane of the circle in
        convenient coordinates for non-tensor linear algebra
        """
        try:

            # if the circle normal is the negative of the xy normal,
            # identity works
            if _np.array_equal(self._normal_vector, _np.array([0, 0, -1])):
                rotation = _np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
                self._rt_to_xy = _rt_to_homog(
                    rotation, _np.transpose(translation)
                )
                translation = -1 * (self._circle_origin)
                self._rt_to_xy = _to_homog(rotation, translation)
                return self._rt_to_xy
            else:
                xy_unit_normal = _np.array([0, 0, 1])
                ortho = _np.cross(self._normal_vector, xy_unit_normal)
                ortho_x = _np.array(
                    [
                        [0, -ortho[2], ortho[1]],
                        [ortho[2], 0, -ortho[0]],
                        [-ortho[1], ortho[0], 0],
                    ]
                )
                rotation = (
                    _np.identity(3)
                    + ortho_x
                    + ortho_x
                    @ ortho_x
                    * (1 / (1 + _np.dot(self._normal_vector, xy_unit_normal)))
                )
                translation = rotation @ (-1 * (self._circle_origin))
                self._rt_to_xy = _to_homog(rotation, translation)
                return self._rt_to_xy
        except AttributeError as ae:
            self.compute_circle()
            self.compute_rt_to_xy()

    def get_radius(self):
        """
        Returns the radius of the circle of allowed atom positions

        Computes the radius, normal, and circle ori if they not yet computed
        """
        try:
            return self._circle_radius
        except AttributeError as no_rad:
            self.compute_circle()
            return self._circle_radius

    def get_circle_origin(self):
        """
        Returns the origin of the circle of allowed atom positions

        Computes the radius, normal, and circle ori if they not yet computed
        """
        try:
            return self._circle_origin
        except AttributeError as no_rad:
            self.compute_circle()
            return self._circle_origin

    def get_normal_vector(self):
        """
        Returns the unit normal vector of the circle of allowed atom positions

        Computes the radius, normal, and circle ori if they not yet computed
        """
        try:
            return self._normal_vector
        except AttributeError as no_rad:
            self.compute_circle()
            return self._normal_vector

    def get_rt_to_xy_plane(self):
        """
        Returns the RT from orgin/normal vector to 0,0,0/0,0,1 (the xy plane)

        Computes the rt, radius, normal, and circle origin if not yet computed
        """
        try:
            return self._rt_to_xy
        except AttributeError as no_rad:
            self.compute_rt_to_xy()
            return self._rt_to_xy

    def closest_point(self, target):
        """
        Returns a numpy vector of the xyz for the closest point on the CAC
        """
        return _closest(
            target,
            self.get_circle_origin(),
            self.get_normal_vector(),
            self.get_radius(),
        )

    def distance_to_closest(self, target):
        """
        returns the distance to the closest point on the circle to target
        """
        return _np.linalg.norm(target - closest_point(target))

    def chi_to_get_closest(self, target):
        """
        Returns the chi angle that gives the closest point on the circle
        """
        closest = closest_point(target)
        atom_xyz = _np.array(
            [self._atom_xyz.x(), self._atom_xyz.y(), self._atom_xyz.z()]
        )
        input_chi = self._chi
        orig_to_atom = atom_xyz - self.get_circle_origin()
        orig_to_closest = closest - self.get_circle_origin()
        delta_theta = vector_theta(orig_to_closest, orig_to_atom)
        new_chi = input_chi + delta_theta
        return new_chi
