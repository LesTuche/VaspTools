
# General imports
import os
import re
from copy import deepcopy

# ASE and materialsproject imports
from mp_api.client import MPRester
from ase.io import read, write
from ase.build import surface
from ase.build.tools import sort
from ase.constraints import FixAtoms
# Vasptools imports
from vasptools.structureopt import StructureOptimization
from vasptools.vasp_recommended_pp import VASP_RECOMMENDED_PP
from drm_copt_estefania_tools import incar_tags_slab, incar_tags_bulk
from typing import Tuple


def currate_contcar(filepath: str):

    # This regex matches “/” followed by one or more decimal digits.
    remove_num = re.compile(r'/[0-9A-Fa-f]+')

    # Read all lines, stripping out “/12345…” whenever it appears:
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        new_lines = []
        for line in f:
            # If there’s a slash‐number sequence (e.g. “/17880443556”), remove it:
            if '/' in line:
                line = remove_num.sub('', line)
            new_lines.append(line)

     # Overwrite the same file with our modified lines:
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(new_lines)


def make_poscar_from_mpid(mp_id, conventional_unite_cell=True, path: str = "") -> str:
    """
    Void function to make a POSCAR file from a Materials Project ID.
    Create a POSCAR file from a Materials Project ID and save it in the current working directory or a specified path.

    #TODO Maybe it is reasonable to have a DB to save the POSCAR files, so that it is not necessary to download them again and again.

    mp_id : str
        The Materials Project ID of the material (e.g., 'mp-1234').
    name : str
        The name of the output POSCAR file (without extension).
        Full path can be provided. If not provided, it will be saved in the current working directory.
        conventional_unit_cell : bool, optional
        If True, the structure will be converted to a conventional unit cell. Default is True. This results in the format if one downloads the POSCAR file from Materials Project manually.

    Returns

    -------
    str
    path to the saved POSCAR file.


    """
    if not path:
        path = os.getcwd()

    if not os.path.exists(path):
        # Create the directory if it does not exist
        print(f"Directory {path} does not exist. Creating it.")
        os.makedirs(path)

    filepath = os.path.join(path, f"{mp_id}.poscar")

    with MPRester() as mpr:
        try:
            structure = mpr.get_structure_by_material_id(
                mp_id, conventional_unit_cell=True)
        except Exception as e:
            print(f"Error fetching structure for MP ID {mp_id}: {e}")
            return ""

        structure.to(filepath, fmt="poscar")

        return filepath


def make_poscars_from_formula(formula: str, path: str = ""):
    """
    Void function to make a POSCAR file from a chemical formula.

    Creates POSCAR file for all materials with the given formula from the Materials Project database, where the structure is expereimentally observed

    Parameters
    ----------

    formula : str
        The chemical formula of the material (e.g., 'Al2O3', 'Pt').

    path : str, optional
        The path where the POSCAR files will be saved. If not provided, it defaults to the current working directory.

    """
    if not path:
        path = os.getcwd()
    if not os.path.exists(path):
        # Create the directory if it does not exist
        print(f"Directory {path} does not exist. Creating it.")
        os.makedirs(path)

    with MPRester() as mpr:

        try:
            structures = mpr.materials.summary.search(
                formula=formula, theoretical=False, fields=["material_id", "formula_pretty"])

        except Exception as e:
            print(f"Error fetching structures for formula {formula}: {e}")
            return []
        paths = []
        for struct in structures:
            print(
                f"Downloading POSCAR for {struct.material_id} with formula {struct.formula_pretty}")
            # Create a POSCAR file for each material
            paths.append(make_poscar_from_mpid(mp_id=struct.material_id))
        return paths


def prepare_bulk_structure(material: str, incar_tags_user: dict, kspacing: float = 0.15, folder_path: str = ""):
    """
    Prepare the bulk structure for a given material.

    Parameters
    ----------
    material : str
        The material formula (e.g., 'Al2O3', 'Pt') or a Materials Project ID (e.g., 'mp-1234').
        if a formula is provided, all experimental structures with this formula will be downloaded from the Materials Project database.
    incar_tags : dict
        A dictionary containing the INCAR tags for the VASP calculations.
    kspacing : float, optional
        The k-point spacing for the calculations. Default is 0.15.
    folder_name : str, optional
        The name of the folder where the input files will be saved. If None, it defaults to the current working directory.

    Returns
    -------
    None
    """

    if not folder_path:
        folder_path = os.getcwd()

    # Create directory if it doesn't exist

    os.makedirs(folder_path, exist_ok=True)
    filepaths = []
    if material.startswith('mp-'):
        # If a Materials Project ID is provided
        filepaths.append(make_poscar_from_mpid(material, path=folder_path))
    else:
        # If a chemical formula is provided, download all experimental structures with this formula
        filepaths = make_poscars_from_formula(
            formula=material, path=folder_path)

    for i, filepath in enumerate(filepaths):

        if not filepath:
            print(f"Failed to create POSCAR for {material}. Skipping.")
            continue

        print(
            f"Preparing {filepath} files for a bulk calculation of {material} for VASP.")  # Read the structure from the POSCAR file

        # get Id from the filepath
        material_id = os.path.basename(filepath).split('.')[0]
        atoms = read(filepath)
        # to avoid overwriting the original dictionary
        incar_tags = deepcopy(incar_tags_bulk)

        # Add or overwrite the default INCAR tags with the provided ones
        incar_tags.update(incar_tags_user)

        job = StructureOptimization(atoms, incar_tags=incar_tags, kspacing=kspacing,  kspacing_definition='vasp',
                                    potcar_dict=VASP_RECOMMENDED_PP, periodicity='3d',
                                    kpointstype='gamma')
        print(
            f"Writing input files for {material} to {folder_path}/bulk_structure_{i+1}")
        job.write_input_files(folder_name=folder_path +
                              f"/bulk_structure_{material_id}")


def prepare_slab_structure(bulk_path: str, miller_indices: Tuple[int, int, int], layers: int = 4, vacuum: int = 8, incar_tags_user: dict = None):
    """
    PRE:

    Expects that bulk_path has been created with prepare_bulk_structure function.
    Expects that a contcar file is present in the bulk_path directory.
    Expects that the bulk_path directory is already created.



    Prepare the slab structure and write VASP input files for a given bulk material.
    Parameters
    ----------
    bulk_path : str
        The path to the bulk structure files (e.g., 'Al2O3_bulk/bulk_structure_3').
    layers : int, optional

        The number of layers in the slab. Default is 4.
    vacuum : int, optional
        The vacuum thickness in Angstroms. Default is 8.
    incar_tags_user : dict, optional
        A dictionary containing the INCAR tags for the VASP calculations. If None, default tags will be used.
    Returns
    -------
    None
    """

    if not os.path.exists(bulk_path + '/CONTCAR'):
        raise FileNotFoundError(f"CONTCAR file  not found in: {bulk_path}")

    # Curate Contcar file
    currate_contcar(bulk_path + '/CONTCAR')

    name = next((name for name in bulk_path.split(
        '/') if "_bulk" in name), None)

    bulk = read(bulk_path + '/CONTCAR')
    slab = surface(bulk, miller_indices, layers=layers,
                   vacuum=vacuum, periodic=True)
    slab = slab.repeat((2, 2, 1))

    slab.translate((0.5, 0.5, 0.0))
    slab.pbc = (True, True, False)
    slab.wrap()

    # Constrain the bottom two layers
    fix_constraint = FixAtoms(
        indices=[atom.index for atom in slab if atom.position[2] < 0.5 * slab.cell[2, 2]])
    slab.set_constraint(fix_constraint)
    incar_tags = deepcopy(incar_tags_slab)
    if incar_tags_user is not None:
        incar_tags.update(incar_tags_user)
    slab = sort(slab)

    job = StructureOptimization(
        atoms=slab,
        incar_tags=incar_tags,
        kspacing=0.15,
        kspacing_definition='vasp',
        kpointstype='gamma',
        potcar_dict=VASP_RECOMMENDED_PP,
        periodicity='2d',
        # Example for Co, adjust as needed
    )

    job.write_input_files(folder_name=bulk_path + f"/{name}_111_slab")


if __name__ == "__main__":
    # Example usage
    material = 'Pt'  # or 'mp-1234' for a specific Materials Project ID
    incar_tags_user = {

        'ISPIN': 2,

    }
    # prepare_bulk_structure(material, incar_tags_user,
    #                      kspacing=0.15, folder_path='Pt_bulk')
    # print(f"Bulk structure preparation for {material} completed.")

    # currate_contcar('Al2O3_bulk/bulk_structure_3/CONTCAR')

    prepare_slab_structure(bulk_path='Pt_bulk/bulk_structure_1',
                           miller_indices=(1, 1, 1), layers=4, vacuum=8,
                           )
