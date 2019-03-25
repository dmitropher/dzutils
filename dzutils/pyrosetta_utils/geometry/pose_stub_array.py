import numpy as _np
from .geometry.homog import stub_to_homog as _stub_to_homog
from .geometry.rt_utils import stub_from_residue as _stub_from_residue


class PoseStubArray(_np.ndarray):
    def __new__(cls, pose=None, seqpos=None, atoms=("CA", "N", "CA", "C")):
        obj = np.asarray(
            _stub_to_homog(_stub_from_residue(pose.residue(seqpos), *atoms))
        ).view(cls)
        obj.pose = pose
        obj.seqpos = seqpos
        obj.atoms = atoms
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.pose = getattr(obj, "pose", None)
        self.seqpos = getattr(obj, "seqpos", None)
        self.atoms = getattr(obj, "atoms", None)
