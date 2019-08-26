def grids_with_aligned_bb(pose, resnum):
    """
    Returns the results of a search to match ca of resnum to a parametric helix

    Works best if its in the neighborhood of "reasonable" helical bundles
    """

    bb_coords = get_bb_coords(pose, resnum)
