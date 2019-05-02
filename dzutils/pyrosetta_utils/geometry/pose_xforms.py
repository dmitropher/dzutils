import numpy as _np
from dzutils.pyrosetta_utils.geometry.homog import (
    stub_to_homog as _stub_to_homog,
)
from dzutils.pyrosetta_utils.geometry.homog import (
    homog_relative_transform as _stubs_to_rt,
)
from dzutils.pyrosetta_utils.geometry.rt_utils import (
    stub_from_residue as _stub_from_residue,
)
import pyrosetta.rosetta as _pyr

# import dzutils.func_utils.MultiMethod as MultiMethod
# import dzutils.func_utils.multimethod as multimethod


class PoseRTArray(_np.ndarray):
    """
    A class designed to wrap a numpy array with pose and residue options

    specifically used to save the relative transform between two residues,
    keeping track of the important atoms.

    Instantiating this object with just stub xforms overrides other behavior,
    otherwise the stubs are computed from the input pose and seqpos
    """

    def __new__(
        cls,
        pose_start=None,
        pose_end=None,
        seqpos_start=None,
        seqpos_end=None,
        stub_xform1=None,
        stub_xform2=None,
        atoms_start=("CA", "N", "CA", "C"),
        atoms_end=("CA", "N", "CA", "C"),
    ):
        if not isinstance(stub_xform1, type(None)) and not isinstance(
            stub_xform2, type(None)
        ):
            obj = _np.asarray(_stubs_to_rt(stub_xform1, stub_xform2)).view(cls)
        else:
            stub_xform1 = _stub_to_homog(
                _stub_from_residue(
                    pose_start.residue(seqpos_start), *atoms_start
                )
            )
            stub_xform2 = _stub_to_homog(
                _stub_from_residue(pose_end.residue(seqpos_end), *atoms_end)
            )
            obj = _np.asarray(_stubs_to_rt(stub_xform1, stub_xform2)).view(cls)
        obj.pose_start = pose_start
        obj.pose_end = pose_end
        obj.seqpos_start = seqpos_start
        obj.seqpos_end = seqpos_end
        obj._atoms_start = atoms_start
        obj._atoms_end = atoms_end
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.pose_start = getattr(obj, "pose_start", None)
        self.pose_end = getattr(obj, "pose_end", None)
        self.seqpos_start = getattr(obj, "seqpos_start", None)
        self.seqpos_end = getattr(obj, "seqpos_end", None)
        self._atoms_start = getattr(obj, "atoms_start", None)
        self._atoms_end = getattr(obj, "atoms_end", None)

    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        pickled_state = super(PoseRTArray, self).__reduce__()
        # Create our own tuple to pass to __setstate__
        new_state = pickled_state[2] + (
            self.pose_start,
            self.pose_end,
            self.seqpos_start,
            self.seqpos_end,
            self._atoms_start,
            self._atoms_end,
        )
        return (pickled_state[0], pickled_state[1], new_state)

    def __setstate__(self, state):
        self._atoms_end = state[-1]  # Set our attributes
        self._atoms_start = state[-2]
        self.seqpos_end = state[-3]
        self.seqpos_start = state[-4]
        self.pose_end = state[-5]
        self.pose_start = state[-6]
        # Call the parent's __setstate__ with the other tuple elements.
        super(PoseRTArray, self).__setstate__(state[0:-6])

    def minimal_pose(self):
        """
        Changes the pose_start and pose_end to a single pose with two residues

        Residue 1 is the residue at seqpos_start 2 is the one at seqpos_end
        """
        pose = _pyr.core.pose.Pose()
        pose.append_residue_by_bond(self.pose_start.residue(self.seqpos_start))
        pose.append_residue_by_jump(self.pose_end.residue(self.seqpos_end), 1)
        self.pose_start = pose
        self.pose_end = pose
        self.seqpos_start = 1
        self.seqpos_end = 2


class PoseStubArray(_np.ndarray):
    """
    A class designed to wrap a numpy array with pose and residue options
    """

    def __new__(cls, pose=None, seqpos=None, atoms=("CA", "N", "CA", "C")):
        obj = _np.asarray(
            _stub_to_homog(_stub_from_residue(pose.residue(seqpos), *atoms))
        ).view(cls)
        obj.pose = pose
        obj.seqpos = seqpos
        obj._atoms = atoms
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.pose = getattr(obj, "pose", None)
        self.seqpos = getattr(obj, "seqpos", None)
        self._atoms = getattr(obj, "atoms", None)

    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        pickled_state = super(PoseStubArray, self).__reduce__()
        # Create our own tuple to pass to __setstate__
        new_state = pickled_state[2] + (self.pose, self.seqpos, self._atoms)
        return (pickled_state[0], pickled_state[1], new_state)

    def __setstate__(self, state):
        self._atoms = state[-1]  # Set our attributes
        self.seqpos = state[-2]
        self.pose = state[-3]
        # Call the parent's __setstate__ with the other tuple elements.
        super(PoseStubArray, self).__setstate__(state[0:-3])

    def set_atoms(self, atoms):
        """
        Takes a tuple of atoms and makes them the atoms for this object.

        Recomputes the h transform
        """
        self._atoms = atoms
        _np.place(
            self,
            _np.ones_like(self),
            _np.asarray(
                _stub_to_homog(
                    _stub_from_residue(self.pose.residue(self.seqpos), *atoms)
                )
            ),
        )

    def get_rt_to_stub(self, stub):
        """
        Returns a wrapper for a homogxform ndarray describing the rt to the stub

        Creates a local rotation translation homogeneous transform
        returns a wrapper object that saves the pose, seqpos, and atoms that the
        rt was generated from.
        """
        rt_array = PoseRTArray(
            pose_start=self.pose,
            pose_end=stub.pose,
            seqpos_start=self.seqpos,
            seqpos_end=stub.seqpos,
            stub_xform1=self,
            stub_xform2=stub,
            atoms_start=self._atoms,
            atoms_end=stub._atoms,
        )
        return rt_array


class RotamerRTArray(_np.ndarray):
    """
    Contains an RT from one set of atoms to another within a single residue

    Also caches the chi, phi, psi. Uses the nomenclature target and base,
    where base is the coordinate frame of the transform and target is
    the R group atoms of interest.


    Defaults to the N,CA,C for the base atoms
    """

    def __new__(
        cls, residue=None, base_atoms=("N", "CA", "C"), target_atoms=None
    ):
        # make the rt from base to target atoms, store it as the 2d array
        obj = _np.asarray(
            _stubs_to_rt(
                _stub_to_homog(_stub_from_residue(residue, *base_atoms)),
                _stub_to_homog(_stub_from_residue(residue, *target_atoms)),
            )
        ).view(cls)
        # generate a blank pose and add our residue into it
        # Residue objects don't update xyz with chi, this is sort of a hack to
        # avoid reimplementing internal coordinate stuff
        pose = _pyr.core.pose.Pose()
        pose.append_residue_by_bond(residue)
        res = pose.residue(1)
        obj.residue = res
        obj._target_atoms = target_atoms
        obj._base_atoms = base_atoms
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.residue = getattr(obj, "residue", None)
        self._target_atoms = getattr(obj, "_target_atoms", None)
        self._base_atoms = getattr(obj, "_base_atoms", None)

    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        pickled_state = super(RotamerRTArray, self).__reduce__()
        # # Create our own tuple to pass to __setstate__
        new_state = pickled_state[2] + (
            self.get_pose,
            self._target_atoms,
            self._base_atoms,
        )
        return (pickled_state[0], pickled_state[1], new_state)

    def __setstate__(self, state):
        self._base_atoms = state[-1]  # Set our attributes
        self._target_atoms = state[-2]
        self.residue = state[-3].residue(1)
        # Call the parent's __setstate__ with the other tuple elements.
        super(RotamerRTArray, self).__setstate__(state[0:-3])

    def _recompute_xform(self):
        """
        Utility function for when the chis are reset
        """
        _np.place(
            self,
            _np.ones_like(self),
            _np.asarray(
                _stubs_to_rt(
                    _stub_to_homog(
                        _stub_from_residue(self.residue, *self._base_atoms)
                    ),
                    _stub_to_homog(
                        _stub_from_residue(self.residue, *self._target_atoms)
                    ),
                )
            ),
        )

    def set_target_atoms(self, atoms):
        """
        Takes a tuple of atoms and makes them the target atoms for this object.

        Recomputes the h transform
        """
        self._target_atoms = atoms
        self._recompute_xform()

    def set_base_atoms(self, atoms):
        """
        Takes a tuple of atoms and makes them the base atoms for this object.

        Recomputes the h transform
        """
        self._base_atoms = atoms
        self._recompute_xform()

    def set_chi(self, chi_num, value):
        """
        Sets the given chi, recomputes xform
        """
        pose = self.get_pose()
        pose.set_chi(chi_num, 1, value)
        self.residue = pose.residue(1)
        self._recompute_xform()

    def set_all_chi(self, *values):
        """
        Sets all chi, recomputes xform
        """
        pose = self.get_pose()
        for chi_num, value in enumerate(values, 1):
            pose.set_chi(chi_num, 1, value)
        self.residue = pose.residue(1)
        self._recompute_xform()

    def get_pose(self):
        """
        Returns a pose where the only residue is the chosen residue
        """
        pose = _pyr.core.pose.Pose()
        pose.append_residue_by_bond(self.residue)
        return pose

    def get_chi(self, chi_num):
        """
        Returns the given chi
        """
        return self.get_pose().chi(chi_num, 1)

    def get_chi_list(self):
        """
        Returns a list of chis in order of chi_num

        Sorry, 0-indexed
        """
        return [
            self.get_chi(chi_num)
            for chi_num in range(1, self.residue.nchi() + 1)
        ]

    def get_pose_stub_array(self, base=True):
        """
        Returns the PoseStubArrray for base by default, target otherwise
        """
        return PoseStubArray(
            pose=self.get_pose(),
            seqpos=1,
            atoms=self._base_atoms if base else self._target_atoms,
        )

    def get_rt_from_base(self, stub_array):
        """
        Returns PoseRTArray ndarray describing the rt from base atoms to stub

        Supports only PoseStubArray objects

        Creates a local rotation translation homogeneous transform
        returns a wrapper object that saves the pose, seqpos, and atoms that the
        rt was generated from.
        """
        # Stuff here
        return self.get_pose_stub_array().get_rt_to_stub(stub_array)

    def get_rt_from_target(self, stub_array):
        """
        Returns an RT from target atoms to stub

        Supports only PoseStubArray objects
        """
        # Stuff here
        return self.get_pose_stub_array(base=False).get_rt_to_stub(stub_array)


def generate_pose_rt(
    pose_1,
    pose_2,
    resnum_1,
    resnum_2,
    atoms_1=("CA", "N", "CA", "C"),
    atoms_2=("CA", "N", "CA", "C"),
    minimal=False,
):
    """
    Convenience wrapper: returns the PoseRTArray for two poses etc

    Written to shorten syntax a bit; defaults to bb atoms
    """
    pose_rt_array = PoseStubArray(
        pose=pose_1, seqpos=resnum_1, atoms=atoms_1
    ).get_rt_to_stub(
        PoseStubArray(pose=pose_2, seqpos=resnum_2, atoms=atoms_2)
    )
    if minimal:
        pose_rt_array.minimal_pose()
        return pose_rt_array
    else:
        return pose_rt_array


def generate_pose_rt_between_res(
    pose,
    resnum_1,
    resnum_2,
    atoms_1=("CA", "N", "CA", "C"),
    atoms_2=("CA", "N", "CA", "C"),
    minimal=False,
):
    """
    Convenience wrapper: defines RT between different residues of the same pose
    """
    pose_rt_array = generate_pose_rt(
        pose, pose, resnum_1, resnum_2, atoms_1, atoms_2
    )
    if minimal:
        pose_rt_array.minimal_pose()
        return pose_rt_array
    else:
        return pose_rt_array
