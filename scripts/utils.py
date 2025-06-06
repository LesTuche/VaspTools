import os
from mp_api.client import MPRester
from ase.io import read, write
from vasptools.structureopt import StructureOptimization
from vasptools.vasp_recommended_pp import VASP_RECOMMENDED_PP
from drm_copt_estefania_tools import incar_tags_bulk
from copy import deepcopy


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
                formula="Pt", theoretical=False, fields=["material_id"])

        except Exception as e:
            print(f"Error fetching structures for formula {formula}: {e}")
            return []
        paths = []
        for struct in structures:
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

    for filepath in filepaths:
        print(
            f"Preparing {filepath} files for a bulk calculation of {material} for VASP.")

        if not filepath:
            print(f"Failed to create POSCAR for {material}. Skipping.")
            continue

        # Read the structure from the POSCAR file

        atoms = read(filepath)
        # to avoid overwriting the original dictionary
        incar_tags = deepcopy(incar_tags_bulk)

        # Add or overwrite the default INCAR tags with the provided ones
        incar_tags.update(incar_tags_user)

        job = StructureOptimization(atoms, incar_tags=incar_tags, kspacing=kspacing,  kspacing_definition='vasp',
                                    potcar_dict=VASP_RECOMMENDED_PP, periodicity='3d',
                                    kpointstype='gamma')
        print(f"Writing input files for {material} to {folder_path}/bulk")
        job.write_input_files(folder_name=folder_path+"/bulk")


if __name__ == "__main__":
    # Example usage
    material = 'Pt'
    incar_tags_user = {

        'ISPIN': 2,

    }
    prepare_bulk_structure(material, incar_tags_user,
                           kspacing=0.15, folder_path='Pt_bulk')
    print(f"Bulk structure preparation for {material} completed.")
