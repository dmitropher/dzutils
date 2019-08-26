# Script to fit ZCON four helix bundles and append fifth and sixth 'buttressing helices' parametrically
# Basile Wicky 190606

# Import libraries
from Bio.PDB import PDBParser  # for parsing input structure
from math import cos, sin, tan, asin, acos, radians, sqrt, degrees, atan
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize, report_fit  # fot fitting helix
import io  # for printing output to file
from contextlib import redirect_stdout  # for printing output to file
import sys

# INPUT STRUCTURE
name = sys.argv[1]

# Switch to TRUE if want movies of minimazation trajectories
make_movies = False


######################## FUNCTIONS ##########################

# Useful for debugging intermediates
def makePDB(coordinates, name):
    resid = 1
    atomnumb = 2
    with open("%s.pdb" % name, "w") as f:
        for i in coordinates:
            f.write(
                "ATOM{:>7s}  CA  ALA A{:>4s}     {:>7.3f} {:>7.3f} {:>7.3f}  1.00  0.00           C\n".format(
                    str(atomnumb),
                    str(resid),
                    float(i[0]),
                    float(i[1]),
                    float(i[2]),
                )
            )
            resid = resid + 1
            atomnumb = atomnumb + 10
        f.write("END")
    return


# Useful for visualizing fitting
def makePDBmovie(coordinates, name):
    with open("%s.pdb" % name, "w") as f:
        model = 1
        for frame in coordinates:
            resid = 1
            atomnumb = 2
            f.write("MODEL{:>8s}\n".format(str(model)))
            for i in frame:
                f.write(
                    "ATOM{:>7s}  CA  ALA A{:>4s}     {:>7.3f} {:>7.3f} {:>7.3f}  1.00  0.00           C\n".format(
                        str(atomnumb),
                        str(resid),
                        float(i[0]),
                        float(i[1]),
                        float(i[2]),
                    )
                )
                resid = resid + 1
                atomnumb = atomnumb + 10
            f.write("TER\n")
            f.write("ENDMDL\n")
            model = model + 1
        f.write("END")
    return


#### FUNCTIONS FOR PARAMETRIC FITTING ####

# Parametric helix equation (from Huang et al. (2014), described in SI)
# Contains small patch to make it resemble the BundleGridSampler mover (define helix 'origin' at the center)

# ----------MAKE CARTESIAN COORDINATES FOR CA OF RESIDUE t---------------
def cartesian(r0, omega0, omega1, phi0, phi1, delta_z, t):

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


# ----------MAKE CARTESIAN COORDINATES FOR HELIX 'TRACE'---------------
def xyz_trace(r0, omega0, phi0, delta_z, t):

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

    # Only first part of equation, corresponds to the middle 'trace' of each helix, without helical rotation
    x = r0 * cos(omega0 * t + phi_prime0)
    y = r0 * sin(omega0 * t + phi_prime0)
    z = ((omega0 * r0) / (tan(alpha))) * t

    return [x, y, z]


# ---------MAKE ARRAY OF XYZ COORDINATES FOR ALL CA-----------------------
def moving(r0, omega0, omega1, phi0, phi1, delta_z, helix_length, invert):

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


# -----------OBJECTIVE RMSD FUNCTION USED DURING MINIMIZATION---------------


def rmsd_array(params, helix_length, target):
    # read parameter dictionary
    r0 = params["r0"]
    omega0 = params["omega0"]
    omega1 = params["omega1"]
    phi0 = params["phi0"]
    phi1 = params["phi1"]
    delta_z = params["delta_z"]
    invert = params["invert"]

    # array of pairwise vectors
    subtract_coord = (
        moving(r0, omega0, omega1, phi0, phi1, delta_z, helix_length, invert)
        - target
    )

    # array of pairwise distances
    rmsd_array = np.sqrt(np.sum(np.power(subtract_coord, 2), axis=1))

    # global rmsd
    rmsd = np.sqrt(
        (1 / helix_length)
        * np.sum(np.sum(np.power(subtract_coord, 2), axis=1))
    )

    # Append coordinates and RMSD of minimization trajectories for visualization
    if chain == "chA":
        rmsd_traj_chA.append(rmsd)
        if make_movies == True:
            movieA.append(
                moving(
                    r0,
                    omega0,
                    omega1,
                    phi0,
                    phi1,
                    delta_z,
                    helix_length,
                    invert,
                )
            )
    if chain == "chB":
        rmsd_traj_chB.append(rmsd)
        if make_movies == True:
            movieB.append(
                moving(
                    r0,
                    omega0,
                    omega1,
                    phi0,
                    phi1,
                    delta_z,
                    helix_length,
                    invert,
                )
            )
    if chain == "chC":
        rmsd_traj_chC.append(rmsd)
        if make_movies == True:
            movieC.append(
                moving(
                    r0,
                    omega0,
                    omega1,
                    phi0,
                    phi1,
                    delta_z,
                    helix_length,
                    invert,
                )
            )
    if chain == "chD":
        rmsd_traj_chD.append(rmsd)
        if make_movies == True:
            movieD.append(
                moving(
                    r0,
                    omega0,
                    omega1,
                    phi0,
                    phi1,
                    delta_z,
                    helix_length,
                    invert,
                )
            )

    return rmsd_array


############ FUNCTIONS FOR GEOMETRICAL TRANSFORMATIONS ###############

# Define translation matrix
def translation(translation_axis):  # general along any arbitrary vector
    return np.array(
        [
            [1, 0, 0, translation_axis[0]],
            [0, 1, 0, translation_axis[1]],
            [0, 0, 1, translation_axis[2]],
            [0, 0, 0, 1],
        ]
    )


# Define translation matrix with distance scaling
def translation_d(
    translation_axis, d
):  # general along any arbitrary vector AND scaled by distance d
    return np.array(
        [
            [1, 0, 0, translation_axis[0] * d],
            [0, 1, 0, translation_axis[1] * d],
            [0, 0, 1, translation_axis[2]],
            [0, 0, 0, 1],
        ]
    )


# Define general rotation matrix (for given angle (in radians)) around arbitrary axis
def rotation(theta, axis):
    # ux^2+uy^2+uz^2 = 1 is required !!!
    ux = axis[0]
    uy = axis[1]
    uz = axis[2]
    return np.array(
        [
            [
                cos(theta) + (ux ** 2) * (1 - cos(theta)),
                ux * uy * (1 - cos(theta)) - uz * sin(theta),
                ux * uz * (1 - cos(theta)) + uy * sin(theta),
            ],
            [
                uy * ux * (1 - cos(theta)) + uz * sin(theta),
                cos(theta) + (uy ** 2) * (1 - cos(theta)),
                uy * uz * (1 - cos(theta)) - ux * sin(theta),
            ],
            [
                uz * ux * (1 - cos(theta)) - uy * sin(theta),
                uz * uy * (1 - cos(theta)) + ux * sin(theta),
                cos(theta) + (uz ** 2) * (1 - cos(theta)),
            ],
        ]
    )


################ FUNCTIONS FOR GENERATION HELIX BACKBONE FROM CA TRACE ################
# From PEILONG with modifications
# convert CA trace to full backbone model by superimposing on ideal template (ideal.pdb) by matching 3 consecutive CA atoms

# Define 'stub' (from 3 atoms) - necessary for computing the transformation matrix that matches the structures
def stub(
    b, a, c
):  # a,b,c are the vectors of a, b, c with respect to the general coordinate frame
    e1 = (a - b) / np.linalg.norm(a - b)
    e3 = np.cross(e1, (c - b)) / np.linalg.norm(np.cross(e1, (c - b)))
    e2 = np.cross(e1, e3) / np.linalg.norm(np.cross(e1, e3))
    partial_matrix = np.array([e1, e2, e3, b]).T
    extra_line = np.array([[0, 0, 0, 1]])
    stub_matrix = np.append(partial_matrix, extra_line, axis=0)

    return stub_matrix


# Read original PDB in order to append the new helix to it
zcon = open(name, "r").readlines()

# Function that takes CA coordinates as input and return ideal helix backbone appended to fitted PDB in one PDB file
def make_BB_from_CA(CA_coordinates, suffix, file_name, buttress):

    # Generate 'ideal' stub
    stub_file = map(
        str.split,
        open(
            "/home/bwicky/Design/add_buttressing_helix_to_zcon/ideal.pdb", "r"
        ).readlines(),
    )

    atom = []
    for line in stub_file:
        atom.append(np.array([float(line[6]), float(line[7]), float(line[8])]))

    ideal_stub = stub(
        atom[6], atom[1], atom[11]
    )  # generate ideal stub from CA coordinates

    with open(
        "buttress_%s/input_structures/%s_%s_%s.pdb"
        % (buttress, buttress, file_name[:-4], suffix),
        "w",
    ) as full_pdb:

        # Write original structure to file
        for line in zcon:
            if line != "END\n":
                full_pdb.write(line)

        # Generate and write bb of generated helix based on its CA trace
        atom_num = 1
        res_num = 0
        chain = "E"
        CA_chain = CA_coordinates[:, 0:3]
        for res in range(1, len(CA_chain) - 1):

            res_num = res_num + 1
            actual_stub = stub(
                CA_chain[res], CA_chain[res - 1], CA_chain[res + 1]
            )  # stub based on CA trace
            transform = np.matmul(
                actual_stub, np.linalg.inv(ideal_stub)
            )  # find transformation matrix between ideal and actual stub

            # N
            coords = np.matmul(transform, np.append(atom[5], 1))
            full_pdb.write(
                "ATOM %6d  N   ALA %s %3d    %8.3f%8.3f%8.3f  1.00  0.00           N\n"
                % (atom_num, chain, res_num, coords[0], coords[1], coords[2])
            )
            atom_num = atom_num + 1

            # CA (use actual CA from trace rather than superimposed one)
            coords = CA_chain[res]
            # If want to use CA from ideal stub instead of actual CA
            # coords=np.matmul(transform,np.append(atom[6],1))
            full_pdb.write(
                "ATOM %6d  CA  ALA %s %3d    %8.3f%8.3f%8.3f  1.00  0.00           C\n"
                % (atom_num, chain, res_num, coords[0], coords[1], coords[2])
            )
            atom_num = atom_num + 1

            # (N)H
            coords = np.matmul(transform, np.append(atom[7], 1))
            full_pdb.write(
                "ATOM %6d  H   ALA %s %3d    %8.3f%8.3f%8.3f  1.00  0.00           H\n"
                % (atom_num, chain, res_num, coords[0], coords[1], coords[2])
            )
            atom_num = atom_num + 1

            # C(O)
            coords = np.dot(transform, np.append(atom[8], 1))
            full_pdb.write(
                "ATOM %6d  C   ALA %s %3d    %8.3f%8.3f%8.3f  1.00  0.00           C\n"
                % (atom_num, chain, res_num, coords[0], coords[1], coords[2])
            )
            atom_num = atom_num + 1

            # O
            coords = np.dot(transform, np.append(atom[9], 1))
            full_pdb.write(
                "ATOM %6d  O   ALA %s %3d    %8.3f%8.3f%8.3f  1.00  0.00           O\n"
                % (atom_num, chain, res_num, coords[0], coords[1], coords[2])
            )
            atom_num = atom_num + 1

        full_pdb.write("END")

    return


def topology_check(
    chA_orientation, chB_orientation, chE_orientation
):  # used for generating correct loop connection variable to pass to the XML
    if (
        chA_orientation == False
        and chB_orientation == True
        and chE_orientation == True
    ):
        variables = ["[E+A,B,C,D]", "[A+B,C,D]"]
    if (
        chA_orientation == False
        and chB_orientation == True
        and chE_orientation == False
    ):
        variables = ["[A,E+B,C,D]", "[B+A,C,D]"]
    if (
        chA_orientation == True
        and chB_orientation == False
        and chE_orientation == True
    ):
        variables = ["[A,E+B,C,D]", "[B+A,C,D]"]
    if (
        chA_orientation == True
        and chB_orientation == False
        and chE_orientation == False
    ):
        variables = ["[E+A,B,C,D]", "[A+B,C,D]"]

    return variables


################## FUNCTIONS END #######################

######### PROTOCOL #############
# Parse input PDB and extract corrdinates of CA
parser = PDBParser()
structure = parser.get_structure("target", name)

chain_lengths = []
ca_target_coordinates = []

num_res = 0
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name() == "CA":
                    ca_target_coordinates.append(atom.get_coord())
            num_res = num_res + 1
        chain_lengths.append(num_res)

targets = {}
targets["chA"] = np.array(ca_target_coordinates[0 : chain_lengths[0]])
targets["chB"] = np.array(
    ca_target_coordinates[chain_lengths[0] : chain_lengths[1]]
)
targets["chC"] = np.array(
    ca_target_coordinates[chain_lengths[1] : chain_lengths[2]]
)
targets["chD"] = np.array(
    ca_target_coordinates[chain_lengths[2] : chain_lengths[3]]
)


# Fit ideal parametric helices to input structure and return corresponding Crick parameters
# ---------- PARAMETERS -------------
# NOT VARIED BY ZIBO
# phi0_guess=radians(180) # FIXED? superhelical phase, i.e. 0, 90, 180, 270 for 4 evenly spaced helices, [degrees] -- BundleGridSampler=delta_omega0
omega0_guess = radians(
    -2.85
)  # FIXED? superhelical twist (-2.85 degrees for two layers), relates to omega1, [degrees] -- BundleGridSampler=omega0
omega1_guess = radians(
    +102.85
)  # FIXED? helical twist (+102.85 degrees for two layers), relates to omega0, [degrees] -- BundleGridSampler=omega1

# VARIED BY ZIBO
r0_guess = 2  # VARY, superhelical radius, [angstrom] -- BundleGridSampler=r0
delta_z_guess = (
    0
)  # VARY, offest along the z axis, [angstrom] -- BundleGridSampler=z0_offset
phi1_guess = radians(
    0
)  # VARY, helical phase (around the internal axis of that helix), [degrees] -- BundleGridSampler=delta_omega1

# For movies of minimization trajectories
movieA = []
movieB = []
movieC = []
movieD = []

# For plotting minimization trajectories
rmsd_traj_chA = []
rmsd_traj_chB = []
rmsd_traj_chC = []
rmsd_traj_chD = []

helix_length_dict = (
    {}
)  # store helix lengths - useful when different helices have very different lengths
fits = {}  # dictonnary that will hold the results from the fits
for chain in targets:

    target_helix = targets[chain]

    # VARIED BY AJASJA
    helix_length = len(target_helix)  # number of residues in the helix, [aa]
    helix_length_dict[chain] = helix_length

    if (target_helix[-1] - target_helix[0])[2] < 0:  # check helix orientation
        invert = True
    else:
        invert = False

    # Generate quick estimates of phi0 to pass as better guess to avoid convergence problems and math domain errors
    if (
        np.average(targets[chain], axis=0)[0] > 0
        and np.average(targets[chain], axis=0)[1] > 0
    ):  # quadrant I
        phi0_guess = radians(
            degrees(
                atan(
                    abs(np.average(targets[chain], axis=0)[1])
                    / abs(np.average(targets[chain], axis=0)[0])
                )
            )
        )
    if (
        np.average(targets[chain], axis=0)[0] < 0
        and np.average(targets[chain], axis=0)[1] > 0
    ):  # quadrant II
        phi0_guess = radians(
            180
            - degrees(
                atan(
                    abs(np.average(targets[chain], axis=0)[1])
                    / abs(np.average(targets[chain], axis=0)[0])
                )
            )
        )
    if (
        np.average(targets[chain], axis=0)[0] < 0
        and np.average(targets[chain], axis=0)[1] < 0
    ):  # quadrant III
        phi0_guess = radians(
            degrees(
                atan(
                    abs(np.average(targets[chain], axis=0)[1])
                    / abs(np.average(targets[chain], axis=0)[0])
                )
            )
            + 180
        )
    if (
        np.average(targets[chain], axis=0)[0] > 0
        and np.average(targets[chain], axis=0)[1] < 0
    ):  # quadrant IV
        phi0_guess = radians(
            360
            - degrees(
                atan(
                    abs(np.average(targets[chain], axis=0)[1])
                    / abs(np.average(targets[chain], axis=0)[0])
                )
            )
        )

    # -----GENERATE PARAMETER DICTIONARY---------
    params = Parameters()
    params.add(
        "r0", value=r0_guess, min=0.001, max=20, vary=True
    )  # avoid negative radii
    params.add(
        "omega0",
        value=omega0_guess,
        min=radians(-4.5),
        max=radians(-0.001),
        vary=True,
    )
    params.add("omega1", value=omega1_guess, vary=True)
    params.add("phi0", value=phi0_guess, vary=True)
    params.add("phi1", value=phi1_guess, vary=True)
    params.add("delta_z", value=delta_z_guess, vary=True)
    params.add("invert", value=invert, vary=False)

    # FIT
    fit = minimize(
        rmsd_array, params, args=(helix_length, target_helix)
    )  # lmfit function
    fits[chain] = fit

# Plot minimization trajectories
plt.plot(rmsd_traj_chA, label="chain A")
plt.plot(rmsd_traj_chB, label="chain B")
plt.plot(rmsd_traj_chC, label="chain C")
plt.plot(rmsd_traj_chD, label="chain D")
plt.title("Minimization trajectories")
plt.xlabel("Minimization steps")
plt.ylabel("RMSD ($\AA$)")
plt.legend(loc="best")
plt.savefig("min_traj.pdf")

# Print full fitting results to file
f = io.StringIO()
with redirect_stdout(f):
    print("-------- CHAIN A ---------")
    report_fit(fits["chA"])

    print("\n\n-------- CHAIN B ---------")
    report_fit(fits["chB"])

    print("\n\n-------- CHAIN C ---------")
    report_fit(fits["chC"])

    print("\n\n-------- CHAIN D ---------")
    report_fit(fits["chD"])
out = f.getvalue()

with open("fit_log.txt", "w") as f_out:
    f_out.write(out)

# Generate file of fit summary
fitA = fits["chA"].params
fitB = fits["chB"].params
fitC = fits["chC"].params
fitD = fits["chD"].params
with open("fit_summary.txt", "w") as f_out:
    f_out.write("\n                 chain A  chain B  chain C  chain D")
    f_out.write(
        "\nRMSD (A)      ={0:>9.3f}{1:>9.3f}{2:>9.3f}{3:>9.3f}".format(
            rmsd_traj_chA[-1],
            rmsd_traj_chB[-1],
            rmsd_traj_chC[-1],
            rmsd_traj_chD[-1],
        )
    )
    f_out.write(
        "\ninvert (Bool) ={0:>9}{1:>9}{2:>9}{3:>9}".format(
            fitA.valuesdict()["invert"],
            fitB.valuesdict()["invert"],
            fitC.valuesdict()["invert"],
            fitD.valuesdict()["invert"],
        )
    )
    f_out.write(
        "\nr0 (A)        ={0:>9.2f}{1:>9.2f}{2:>9.2f}{3:>9.2f}".format(
            fitA.valuesdict()["r0"],
            fitB.valuesdict()["r0"],
            fitC.valuesdict()["r0"],
            fitD.valuesdict()["r0"],
        )
    )
    f_out.write(
        "\ndelta_z (A)   ={0:>9.2f}{1:>9.2f}{2:>9.2f}{3:>9.2f}".format(
            fitA.valuesdict()["delta_z"],
            fitB.valuesdict()["delta_z"],
            fitC.valuesdict()["delta_z"],
            fitD.valuesdict()["delta_z"],
        )
    )
    f_out.write(
        "\nomega0 (deg)  ={0:>9.2f}{1:>9.2f}{2:>9.2f}{3:>9.2f}".format(
            degrees(fitA.valuesdict()["omega0"]),
            degrees(fitB.valuesdict()["omega0"]),
            degrees(fitC.valuesdict()["omega0"]),
            degrees(fitD.valuesdict()["omega0"]),
        )
    )
    f_out.write(
        "\nomega1 (deg)  ={0:>9.2f}{1:>9.2f}{2:>9.2f}{3:>9.2f}".format(
            degrees(fitA.valuesdict()["omega1"]),
            degrees(fitB.valuesdict()["omega1"]),
            degrees(fitC.valuesdict()["omega1"]),
            degrees(fitD.valuesdict()["omega1"]),
        )
    )
    f_out.write(
        "\nphi0 (deg)    ={0:>9.2f}{1:>9.2f}{2:>9.2f}{3:>9.2f}".format(
            degrees(fitA.valuesdict()["phi0"]),
            degrees(fitB.valuesdict()["phi0"]),
            degrees(fitC.valuesdict()["phi0"]),
            degrees(fitD.valuesdict()["phi0"]),
        )
    )
    f_out.write(
        "\nphi1 (deg)    ={0:>9.2f}{1:>9.2f}{2:>9.2f}{3:>9.2f}\n".format(
            degrees(fitA.valuesdict()["phi1"]),
            degrees(fitB.valuesdict()["phi1"]),
            degrees(fitC.valuesdict()["phi1"]),
            degrees(fitD.valuesdict()["phi1"]),
        )
    )

# Generate movie PDB of fitting trajectories
if make_movies == True:
    makePDBmovie(movieA, "movieA")
    makePDBmovie(movieB, "movieB")
    makePDBmovie(movieC, "movieC")
    makePDBmovie(movieD, "movieD")


########### PART WITH GEOMETRICAL TRANSFORMATIONS TO BUTTRESS NEW HELIX TO EXISTING BUNDLE #############

# Generate grid-sampled helices buttressing A and B as defined by user
def generator(
    r0_A,
    r0_B,
    r0_C,
    phi0_A,
    phi0_B,
    phi0_C,
    omega0_A,
    omega0_B,
    omega1_A,
    omega1_B,
    chA_orientation,
    chB_orientation,
    delta_z_A,
    delta_z_B,
    helix_lenght_A,
    helix_length_B,
    buttress,
):

    # Find xy coordinates based on fitted parameteres
    xy_A = np.array([cos(phi0_A) * r0_A, sin(phi0_A) * r0_A])
    xy_B = np.array([cos(phi0_B) * r0_B, sin(phi0_B) * r0_B])
    xy_C = np.array([cos(phi0_C) * r0_C, sin(phi0_C) * r0_C])

    # Compute averages for the parameters that will not sampled on the grid
    omega0_avg = (omega0_A + omega0_B) / 2
    omega1_avg = (omega1_A + omega1_B) / 2
    delta_z_avg = (delta_z_A + delta_z_B) / 2
    helix_length_avg = (helix_lenght_A + helix_length_B) / 2

    # Find the point that is equidistant from A and B on the line that connects the two
    xy_midpoint = ((xy_B - xy_A) / 2) + xy_A
    r_midpoint = np.sqrt(np.sum(np.power(xy_midpoint, 2)))
    if xy_midpoint[0] > 0 and xy_midpoint[1] > 0:  # quadrant I
        angle_midpoint = atan(xy_midpoint[1] / xy_midpoint[0])
    if xy_midpoint[0] < 0 and xy_midpoint[1] > 0:  # quadrant II
        angle_midpoint = atan(xy_midpoint[1] / xy_midpoint[0]) + radians(180)
    if xy_midpoint[0] < 0 and xy_midpoint[1] < 0:  # quadrant III
        angle_midpoint = atan(xy_midpoint[1] / xy_midpoint[0]) + radians(180)
    if xy_midpoint[0] > 0 and xy_midpoint[1] < 0:  # quadrant IV
        angle_midpoint = atan(xy_midpoint[1] / xy_midpoint[0]) + radians(360)

    # ----------- DOT TEST ------------
    # Necessary to ensure that the translation of the newly generated helix is outwards with respect to the bundle
    AvB = np.append(xy_B - xy_A, 0)  # vector AB
    AvC = np.append(
        xy_C - xy_A, 0
    )  # vector AC (helix C or D could be used, ensures that projection is away from it)
    dot_test = False
    if np.cross(AvB, AvC)[2] < 0:
        dot_test = True

    # -------- ARRAYS OF PERPENDICULAR VECTORS AT THE MIDPOINT BETWEEN A AND B ----------
    trace_chA = []
    trace_chB = []
    trace_midpoint = []
    perp_vectors = []
    interation = np.arange(-helix_length_avg / 2, helix_length_avg / 2, 1)
    for k in interation:
        # Compute helical 'trace'
        trace_A = np.array(
            xyz_trace(r0_A, omega0_A, phi0_A, delta_z_avg, k)
        ) + np.array([0, 0, delta_z_avg])
        trace_B = np.array(
            xyz_trace(r0_B, omega0_B, phi0_B, delta_z_avg, k)
        ) + np.array([0, 0, delta_z_avg])

        # Find points that are equidistant from A and B on the line connecting the two
        trace_mid = ((trace_B - trace_A) / 2) + trace_A

        trace_chA.append(trace_A)
        trace_chB.append(trace_B)
        trace_midpoint.append(trace_mid)

        # Find perpendicular vector to the line (with origin at the midpoint) that projects AWAY from the bundle
        if dot_test == False:
            perp = (
                np.cross(trace_B - trace_mid, np.array([0, 0, 1])) + trace_mid
            )
        if dot_test == True:
            perp = (
                np.cross(trace_B - trace_mid, np.array([0, 0, -1])) + trace_mid
            )

        norm_perp_vectors = (perp - trace_mid) / np.linalg.norm(
            perp - trace_mid
        )
        perp_vectors.append(norm_perp_vectors)

        # Useful for visualization
    #         line_v=[trace_A,trace_B] # for visualization of AB lines in PyMOL
    #         makePDB(line_v,'line_%s'%k)
    #         perp_v=[trace_mid,7*norm_perp_vectors+trace_mid] # for visualization of perpendicular vectors in PyMOL
    #         makePDB(perp_v,'perp_%s'%k)

    #     Useful for visualization
    #     makePDB(trace_chA,'trace_chA')
    #     makePDB(trace_chB,'trace_chB')
    #     makePDB(trace_midpoint,'trace_midpoint')

    ################ GRID SAMPLER ##################
    # DOFS: rotation of helix around itself (phi1) / distance from AB (d) / tilting around midle of helix (tilt)
    # ------------ GRID ---------------
    orient_sample = [False, True]
    d_sample = np.arange(8.5, 11.5, 1)
    tilt_sample = np.arange(radians(-20), radians(10), radians(10))
    phi1_sample = np.arange(radians(0), radians(360), radians(40))

    pdb_names = []  # for generating task array
    loop1 = (
        []
    )  # variable to pass to the XML based on the connectivity that is possible for building the first loop
    loop2 = (
        []
    )  # variable to pass to the XML based on the connectivity that is possible for building the second loop
    for orient in orient_sample:
        for d in d_sample:
            for tilting_offset in tilt_sample:
                for phi1 in phi1_sample:
                    # Generate helix using averaged parameters
                    # +2 for helix length necessary for reconstituation of full backbone from CA later on (otherwise will be shorter)
                    # moving(r0, omega0, omega1, phi0, phi1, delta_z, helix_length,invert)
                    generated = moving(
                        r_midpoint,
                        radians(50),
                        radians(50),
                        angle_midpoint,
                        phi1,
                        delta_z_avg,
                        helix_length_avg + 2,
                        orient,
                    )

                    # Translate middle of helix into xy plane before rotation to avoid lever-arm effects
                    closest_to_origin = np.argmin(
                        np.absolute(np.array(trace_midpoint)), axis=0
                    )[2]
                    shift_v = (
                        trace_midpoint[closest_to_origin]
                        - trace_midpoint[int(len(interation) / 2)]
                    )
                    shift_matrix = translation(shift_v)
                    shifted_list = []
                    for line in generated:
                        add = np.append(line, 1)
                        coord = np.dot(shift_matrix, add)
                        shifted_list.append(coord)
                    shifted = np.array(shifted_list)

                    # Rotate around projection vector
                    tilt_v = trace_midpoint[-1] - trace_midpoint[0]
                    tilt_angle = acos(
                        np.dot(tilt_v, np.array([0, 0, 1]))
                        / (np.linalg.norm(tilt_v))
                    )

                    rotation_matrix = rotation(
                        -tilt_angle + tilting_offset,
                        perp_vectors[int(len(interation) / 2)],
                    )
                    rotated = np.dot(shifted[:, 0:3], rotation_matrix)

                    # Shift back into position
                    shift_back_matrix = translation(-shift_v)
                    shifted_back_list = []
                    for line in rotated:
                        add = np.append(line, 1)
                        coord = np.dot(shift_back_matrix, add)
                        shifted_back_list.append(coord)
                    shifted_back = np.array(shifted_back_list)

                    # Translate outwards along projection vector
                    translate_out_matrix = translation_d(
                        perp_vectors[int(len(interation) / 2)], d
                    )
                    trans_out_list = []
                    for line in shifted_back:
                        coord = np.dot(translate_out_matrix, line)
                        trans_out_list.append(coord)
                    trans_out = np.array(trans_out_list)

                    # Keep track of the parameters in the file name
                    suffix = "_%d_%.1f_%d_%d" % (
                        orient,
                        d,
                        degrees(tilting_offset),
                        degrees(phi1),
                    )

                    # Make the full BB form the CA trace
                    make_BB_from_CA(trans_out, suffix, name, buttress)

                    # Keep track of a few things to generate task list
                    pdb_names.append(
                        "input_structures/%s_%s_%s.pdb"
                        % (buttress, name[:-4], suffix)
                    )
                    chE_orientation = orient
                    topology = topology_check(
                        chA_orientation, chB_orientation, chE_orientation
                    )
                    loop1.append(topology[0])
                    loop2.append(topology[1])

                    # Generates PDBs of intermediate transformations - useful for debugging
    #                     makePDB(generated,'generated_%s'%buttress)
    #                     makePDB(shifted,'shifted_%s'%buttress)
    #                     makePDB(rotated,'rotated_%s'%buttress)
    #                     makePDB(shifted_back,'shifted_back_%s'%buttress)
    #                     makePDB(trans_out,'trans_out_%s'%buttress)

    # Generate command task file to pass to slurm
    rosetta_path = "/software/rosetta/versions/v2019.21-dev60746/bin/rosetta_scripts.hdf5.linuxgccrelease"  # static version
    commandList = []
    for i in range(len(pdb_names)):
        command = (
            "%s @flags -in:file:s %s -parser:script_vars loop1_connectivity=%s loop2_connectivity=%s"
            % (rosetta_path, pdb_names[i], loop1[i], loop2[i])
        )
        commandList.append(command)

    with open("buttress_%s/tasks" % buttress, "w") as f:
        for task in commandList:
            f.write(task + "\n")

    return


###### GENERATE BUTTRESS HELIX TO AB
buttress_in = "AB"
# Extract parameters form fit
r0_A_in = fitA.valuesdict()["r0"]
r0_B_in = fitB.valuesdict()["r0"]
r0_C_in = fitC.valuesdict()["r0"]

phi0_A_in = fitA.valuesdict()["phi0"]
phi0_B_in = fitB.valuesdict()["phi0"]
phi0_C_in = fitC.valuesdict()["phi0"]

omega0_A_in = fitA.valuesdict()["omega0"]
omega0_B_in = fitB.valuesdict()["omega0"]

omega1_A_in = fitA.valuesdict()["omega1"]
omega1_B_in = fitB.valuesdict()["omega1"]

chA_orientation_in = fitA.valuesdict()["invert"]
chB_orientation_in = fitB.valuesdict()["invert"]

delta_z_A_in = fitA.valuesdict()["delta_z"]
delta_z_B_in = fitB.valuesdict()["delta_z"]

helix_lenght_A_in = helix_length_dict["chA"]
helix_length_B_in = helix_length_dict["chB"]

# Generate samples of helices buttresing AB
generator(
    r0_A_in,
    r0_B_in,
    r0_C_in,
    phi0_A_in,
    phi0_B_in,
    phi0_C_in,
    omega0_A_in,
    omega0_B_in,
    omega1_A_in,
    omega1_B_in,
    chA_orientation_in,
    chB_orientation_in,
    delta_z_A_in,
    delta_z_B_in,
    helix_lenght_A_in,
    helix_length_B_in,
    buttress_in,
)

###### GENERATE BUTTRESS HELIX TO CD
buttress_in = "CD"
# Extract parameters form fit
r0_A_in = fitC.valuesdict()["r0"]
r0_B_in = fitD.valuesdict()["r0"]
r0_C_in = fitA.valuesdict()["r0"]

phi0_A_in = fitC.valuesdict()["phi0"]
phi0_B_in = fitD.valuesdict()["phi0"]
phi0_C_in = fitA.valuesdict()["phi0"]

omega0_A_in = fitC.valuesdict()["omega0"]
omega0_B_in = fitD.valuesdict()["omega0"]

omega1_A_in = fitC.valuesdict()["omega1"]
omega1_B_in = fitD.valuesdict()["omega1"]

chA_orientation_in = fitC.valuesdict()["invert"]
chB_orientation_in = fitD.valuesdict()["invert"]

delta_z_A_in = fitC.valuesdict()["delta_z"]
delta_z_B_in = fitD.valuesdict()["delta_z"]

helix_lenght_A_in = helix_length_dict["chC"]
helix_length_B_in = helix_length_dict["chD"]

# Generate samples of helices buttressing CD
generator(
    r0_A_in,
    r0_B_in,
    r0_C_in,
    phi0_A_in,
    phi0_B_in,
    phi0_C_in,
    omega0_A_in,
    omega0_B_in,
    omega1_A_in,
    omega1_B_in,
    chA_orientation_in,
    chB_orientation_in,
    delta_z_A_in,
    delta_z_B_in,
    helix_lenght_A_in,
    helix_length_B_in,
    buttress_in,
)
