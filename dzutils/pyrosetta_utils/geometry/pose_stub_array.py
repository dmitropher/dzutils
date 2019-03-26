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
        Returns a homogenous transform ndarray describing the rt to the stub

        Creates a local rotation translation homogeneous transform
        """
        return _stubs_to_rt(self, stub)
