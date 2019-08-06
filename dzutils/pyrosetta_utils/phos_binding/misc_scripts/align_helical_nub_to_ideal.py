from dzutils.pyrosetta_utils.geometry.parametric import ca_array


def align_frag_to_ideal_helix(
    frag, align_stub, r0, omega0, omega1, phi0, phi1, delta_z, invert
):
    """
    Align and return frag by stub to ideal helix described by the params given
    """
    ca_array(r0, omega0, omega1, phi0, phi1, delta_z, 3, invert)


def main():
    """
    """


if __name__ == "__main__":
    main()
