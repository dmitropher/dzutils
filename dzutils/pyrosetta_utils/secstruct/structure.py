# These classes are meant as containers for general information about pieces of
# secondary structure, as well as utility methods for finding and assessing the
# continued chemical properties of these elements.
import pyrosetta.rosetta as _pyr
from itertools import groupby as _gb


class SecondaryStructureContainerFactory:
    def __init__(self):
        self._creators = {}

    def register_format(self, dssp_type, creator):
        self._creators[dssp_type] = creator

    def get_container(self, pose, start_pos, end_pos, dssp_type):
        creator = self._creators.get(dssp_type)
        if not creator:
            return SecondaryStructureResidueContainer(
                pose, start_pos, end_pos, dssp_type
            )
        return creator(pose, start_pos, end_pos, dssp_type)


class SecondaryStructureResidueContainer(object):
    """
    An informational container for rosetta pose elements of secondary structure

    Designed to be subclassed into "helix", "strand","gamma turn", etc etc

    This container is not designed to be dynamic and continue to store imformation about the structure over time! Subclasses should have methods to determine continued accuracy of the container for the pose to which it is bound.
    #TODO apply_label()
    #TODO NotImplementedError() for check_residues
    #TODO create_subpose()
    """

    def __init__(self, pose, start_pos, end_pos, dssp_type):
        self.pose = pose
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.residues = pose.residues[start_pos:end_pos]
        self.dssp_type = dssp_type

    def get_range(self):
        """
        Returns a tuple (start_pos,end_pos)
        """
        return self.start_pos, self.end_pos

    def frame_shift(self, n):
        """
        Add n to both the start and end of the SecondaryStructure (accepts neg)
        """
        self.start_pos, self.end_pos = self.start_pos + n, self.end_pos + n

    # TODO: split_structure method, where the structure can be broken at some residue into two identical structures
    def split_structure(self, index):
        """
        return as tuple two strucures: (start,index) & (index +1,end)

        Should error if you try to split somewhere ~unusual~
        """
        return (
            SecondaryStructureResidueContainer(self.start_pos, index),
            SecondaryStructureResidueContainer(index + 1, self.end_pos),
        )


class HelixResidueContainer(SecondaryStructureResidueContainer):
    """
    A container to describe consecutive helical residues

    TODO: support for kinks, pitch, coiled coil pairs/indexing, check helicity
    """

    def __init__(self, pose, start_pos, end_pos, dssp_type):
        super().__init__(pose, start_pos, end_pos, dssp_type)

    # TODO check_helix
    """
    Checks if all residues start-end are still helical
    """


class LoopResidueContainer(SecondaryStructureResidueContainer):
    """
    A secondary structure container for holding all types of loops

    TODO preceding and following attributes that are bound to secondary structure up and down stream of the loop
    TODO internal container of all structured turns in the loop
    """

    def __init__(self, pose, start_pos, end_pos, dssp_type):
        super().__init__(pose, start_pos, end_pos, dssp_type)


# TODO class Turn(LoopResidueContainer)


class StrandResidueContainer(SecondaryStructureResidueContainer):
    """
    A container to describe consecutive beta strand residues

    TODO: support for bulge, pitch, neighboring strands, check beta
    """

    def __init__(self, pose, start_pos, end_pos, dssp_type):
        super().__init__(pose, start_pos, end_pos, dssp_type)

    # TODO check_beta
    """
    Checks if all residues start-end are still beta character
    """


def get_container_creator():
    factory = SecondaryStructureContainerFactory()
    factory.register_format("E", StrandResidueContainer)
    factory.register_format("H", HelixResidueContainer)
    factory.register_format("L", LoopResidueContainer)
    return factory


def get_allowed_dssp_values():
    """
    Returns the DSSP string values registered with the current creator
    """
    from copy import deepcopy

    return deepcopy(get_container_creator()._creators)


def parse_structure_from_dssp(pose, *dssp_types):
    """
    returns a list of secondary structure containers for the given pose

    leaving dssp_types blank allows all types. Unsupported dssp values use the
    generic SecondaryStructureResidueContainer
    """

    dssp_str = _pyr.core.scoring.dssp.Dssp(pose).get_dssp_secstruct()
    creator = get_container_creator()
    return [
        creator.get_container(pose, run[0], run[-1], str(res_type_string))
        for res_type_string, iterator in _gb(
            enumerate(dssp_str, 1), lambda x: x[1]
        )
        for run in ([n for n, p in iterator],)
        if (not dssp_types) or res_type_string in dssp_types
    ]
