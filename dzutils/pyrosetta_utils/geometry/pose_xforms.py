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
