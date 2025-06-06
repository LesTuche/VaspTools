import os
import sys
import numpy as np
from ase.io import read

#main_path = '/Users/hprats/PycharmProjects/collaborations/drm_copt_estefania'

incar_tags_bulk = {
    # ionic relaxation
    'ISIF': 3,
    'IBRION': 2,
    'NSW': 500,
    # convergence criteria
    'EDIFF': 1E-06,
    'EDIFFG': -0.001,
    'NELM': 300,
    # general electronic settings
    'GGA': 'PE',
    'IVDW': 12,
    'ENCUT': 520,
    'ISPIN': 1,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'LREAL': '.FALSE.',
    'LASPH': '.TRUE.',
    'PREC': 'Accurate',
    # wavefunction and charge density files
    'LWAVE': '.FALSE.',  # WAVECAR
    'LCHARG': '.FALSE.',  # CHGCAR and CHG
    # parallelization
    'NCORE': 16,
    'KPAR': 2,
}


incar_tags_slab = {
    # ionic relaxation
    'IBRION': 2,
    'NSW': 500,
    # convergence criteria
    'EDIFF': 1E-05,
    'EDIFFG': -0.01,
    'NELM': 300,
    # general electronic settings
    'GGA': 'PE',
    'IVDW': 12,
    'ENCUT': 400,
    'ISPIN': 1,
    'ISMEAR': 1,
    'ALGO': 'Fast',
    'LREAL': 'Auto',
    'LASPH': '.TRUE.',
    # dipole corrections
    'LDIPOL': '.TRUE.',
    'IDIPOL': '3',
    'DIPOL': '0.5 0.5 0.5',
    # wavefunction and charge density files
    'LWAVE': '.FALSE.',
    'LCHARG': '.FALSE.',
    # parallelization
    'NCORE': 12,
    'KPAR': 2
}


incar_tags_mlneb = {
    # convergence criteria
    'ediff': '1e-05',
    'ediffg': -0.01,
    'nelm': 300,
    # general electronic settings
    'xc': 'pbe',
    'gga': 'PE',
    'ivdw': 12,
    'encut': 400,
    'ispin': 1,
    'ismear': 1,
    'algo': 'Fast',
    'lreal': 'Auto',
    'lasph': 'True',
    # dipole corrections
    'ldipol': 'True',
    'idipol': 3,
    'dipol': '[0.5, 0.5, 0.5]',
    # wavefunction and charge density files
    'lwave': 'False',
    'lcharg': 'False',
    # parallelization
    'kpar': 2
}


def get_elements(formula):
    num_h, num_o, num_c = 0, 0, 0
    for i in range(len(formula)):
        if formula[i] == 'C':
            num_c += 1
        elif formula[i] == 'O':
            num_o += 1
        elif formula[i] == 'H':
            num_h += 1
        else:
            if formula[i - 1] == 'C':
                num_c += int(formula[i]) - 1
            if formula[i - 1] == 'O':
                num_o += int(formula[i]) - 1
            if formula[i - 1] == 'H':
                num_h += int(formula[i]) - 1
    return num_h, num_o, num_c


def get_formation_energy(formula, total_energy, energy_slab):

    reference_energies = {}
    for gas_molecule in ['H2', 'CH4', 'CO2']:
        atoms = read(f'{main_path}/gas/{gas_molecule}/vasprun.xml', index=-1)
        reference_energies[gas_molecule] = atoms.get_potential_energy()

    ref_h = reference_energies['H2'] / 2
    ref_c = reference_energies['CH4'] - 4 * ref_h
    ref_o = (reference_energies['CO2'] - ref_c) / 2

    num_h, num_o, num_c = get_elements(formula)
    formation_energy = total_energy - energy_slab - num_h * ref_h - num_o * ref_o - num_c * ref_c

    return formation_energy


def get_energy_mlneb(path_mlneb):
    # Check if ibrion1 exists
    path_ibrion1 = f"{path_mlneb}/dimer/ibrion1"
    if os.path.isdir(path_ibrion1):
        return read(f"{path_ibrion1}/vasprun.xml").get_total_energy()

    # If not, check if dimer exists
    path_dimer = f"{path_mlneb}/dimer"
    if os.path.isdir(path_dimer):
        return read(f"{path_dimer}/vasprun.xml").get_total_energy()

    # If not, check if ML-NEB exists, else use last_predicted_path.traj
    for fname in ["ML-NEB.traj", "last_predicted_path.traj"]:
        traj_path = os.path.join(path_mlneb, fname)
        if os.path.isfile(traj_path):
            images = read(traj_path, index=":")
            return max(img.get_potential_energy() for img in images)

    sys.exit(f"{path_mlneb}: Trajectory file for TS not found")


def get_vib_list(path, ignore_modes=0):
    # ignore_modes = 0: all modes are vibrations
    # ignore_modes = 1: transition state for surface processes
    # ignore_modes = 5: gas-phase linear molecule
    # ignore_modes = 6: gas-phase nonlinear molecule

    if os.path.isdir(f"{path}/dimer"):
        vib_path = f"{path}/dimer/vibrations"
    else:
        vib_path = f"{path}/vibrations"
    with open(f"{vib_path}/vibrations.txt") as infile:
        lines = infile.readlines()

    # Only process the expected lines (ignoring the last 'ignore_modes' lines)
    valid_lines = lines[:len(lines) - ignore_modes]

    vib_values = []
    for line in valid_lines:
        if 'f/i' in line:
            print(f"Warning: extra imaginary mode found in {vib_path}/vibrations.txt, ignored")
            continue  # Skip this mode entirely
        vib_values.append(float(line.split()[-2]) / 1000.0)

    return np.array(vib_values)  # in eV


