from copy import deepcopy
from ase.io import read
from ase.lattice.tetragonal import SimpleTetragonalFactory
from vasptools.structureopt import StructureOptimization
from vasptools.vasp_recommended_pp import VASP_RECOMMENDED_PP
from drm_copt_estefania.drm_copt_estefania_tools import  incar_tags_bulk
from utils import make_poscar_from_mp


main_path = os.getcwd()

"""
TODO: - Make a importable function out of this script to prepare the bulk structures.
    : - Input should be a list of materials and their corresponding Materials Project IDs.
    : - Output should be a folder with the input files for each material.

    Fucntion which fetches the CIF of the material from Materials Project. It should be placed in a folder bulk_structures, and a set should be created to check wheter we need to refetch it. 




"""





def prepare_bulk_structure(material, incar_tags, kspacing=0.15, folder_name=None):




    """

    materials: dict
        A dictionary where keys are material symbols (e.g., 'Co', 'Pt') and values are their corresponding Materials Project IDs.
    incar_tags: dict
        A dictionary containing the INCAR tags for the VASP calculations. They will be added or overwirte the default ones.
    kspacing: float
        The k-point spacing for the calculations. Uses the VASP convention.
    folder_name: str, optional
        The name of the folder where the input files will be saved. If None, it defaults to the current working directory.


    
    """
    
    if folder_name is None:
        folder_name=os.getcwd()


    

    

    pass
    



materials = {
    #'Co': 'mp-102',
    'Pt': 'mp-126',
    #  'CoPt': 'mp-949'   # problem: it is the primitive cell, prepare the conventional manually
}

for material in materials:

    if f"main_path/bulk/" not in os.listdir(f'{main_path}'):
        os.mkdir(f'{main_path}/bulk')


    atoms = read(f'{main_path}/bulk/{material}.poscar')

    incar_tags = deepcopy(incar_tags_bulk)

    if 'Co' or "Pt" in material:
        incar_tags['ISPIN'] = 2

    job = StructureOptimization(
        atoms=atoms,
        incar_tags=incar_tags,
        kspacing=0.15,
        kspacing_definition='vasp',
        kpointstype='gamma',
        potcar_dict=VASP_RECOMMENDED_PP,
        periodicity='3d'
    )
    job.write_input_files(folder_name=f'{main_path}/bulk/{material}')


class CoPtFactory(SimpleTetragonalFactory):
    bravais_basis = [[0.0, 0, 0.0], [0.5, 0.5, 0.0], [0.5, 0, 0.5], [0.0, 0.5, 0.5]]
    element_basis = (0, 0, 1, 1)  # 0 for Co, 1 for Pt

atoms = CoPtFactory()
atoms = atoms(symbol=('Co', 'Pt'), latticeconstant={'a': 3.859, 'c': 3.694}, size=(1, 1, 1))

incar_tags = deepcopy(incar_tags_bulk)
incar_tags['ISPIN'] = 2

job = StructureOptimization(
    atoms=atoms,
    incar_tags=incar_tags,
    kspacing=0.66,
    kpointstype='gamma',
    potcar_dict=VASP_RECOMMENDED_PP,
    periodicity='3d'
)
job.write_input_files(folder_name=f'{main_path}/bulk/CoPt')