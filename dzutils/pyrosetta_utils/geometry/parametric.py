# functions derived from Basile Wicky's work on generating parametric buttress helices

from math import cos, sin, tan, asin, acos, radians, sqrt, degrees, atan
import numpy as np


def cartesian(r0, omega0, omega1, phi0, phi1, delta_z, t):
    """
    Returns a list [x,y,z] of cartesian coordinates derived from helical params 
    """

    # SHOULD NOT BE VARIED IF WANT TO RETAIN IDEAL ALPHA-HELIX
    d = (
        1.51
    )  # FIXED, distance between successive residues along the helical axis, [angstrom] -- BundleGridSampler=z1
    r1 = (
        2.26
    )  # FIXED, helical radius, [angstrom] -- BundleGridSampler=r1_peratom

    # ONLY FUNCTIONS OF OTHER PARAMETERS
    alpha = asin(
        (r0 * omega0) / d
    )  # CONSTRAINED (function of other variables), pitch angle, [radians] -- BundleGridSampler=delta_omega1_peratom?
    # N.B. alpha is usually the issue if the script crashes due to a math domain error
    phi_prime0 = (
        phi0 + delta_z * tan(alpha) / r0
    )  # CONSTRAINED (function of other variables), superhelical phase decoupled from delta_z, [radians]

    x = (
        r0 * cos(omega0 * t + phi_prime0)
        + r1 * cos(omega0 * t + phi_prime0) * cos(omega1 * t + phi1)
        - r1
        * cos(alpha)
        * sin(omega0 * t + phi_prime0)
        * sin(omega1 * t + phi1)
    )
    y = (
        r0 * sin(omega0 * t + phi_prime0)
        + r1 * sin(omega0 * t + phi_prime0) * cos(omega1 * t + phi1)
        + r1
        * cos(alpha)
        * cos(omega0 * t + phi_prime0)
        * sin(omega1 * t + phi1)
    )
    z = (
        ((omega0 * r0) / (tan(alpha))) * t
        - r1 * sin(alpha) * sin(omega1 * t + phi1)
        + delta_z
    )

    return [x, y, z]


# ---------MAKE ARRAY OF XYZ COORDINATES FOR ALL CA-----------------------
def ca_array(r0, omega0, omega1, phi0, phi1, delta_z, helix_length, invert):
    """
    Returns an np array of c-alpha coordinates for an ideal helix
    """

    # PATCH TO BRIDGE DIFFERENCES IN HOW THE 'ORIGIN' IS DEFINED
    delta_t = int(
        helix_length / 2
    )  # define an offset of half-helix length (in number of residues) -- BundleGridSampler=delta_t
    # 're-number' indices +/- around middle of helix
    # to patch Vikram's convention (start from middle of helix) and Huang's convention (start at resid 1)

    # Correct for helices that have odd numbers of residues (otherwise fitting helix will be one residue short)
    if (helix_length % 2) == 0:
        residue_renumber_indices = np.arange(-delta_t, +delta_t, 1)
    if (helix_length % 2) != 0:
        residue_renumber_indices = np.arange(-delta_t, +delta_t + 1, 1)

    if invert == True:  # change direction of helix
        residue_renumber_indices = -1 * residue_renumber_indices

    moving_coordinates = []
    for t in residue_renumber_indices:

        xyz = cartesian(r0, omega0, omega1, phi0, phi1, delta_z, t)
        moving_coordinates.append(xyz)

    return np.array(moving_coordinates)
