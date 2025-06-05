from copy import deepcopy
from ase.io import read
from ase.lattice.tetragonal import SimpleTetragonalFactory
from vasptools.structureopt import StructureOptimization
from vasptools.vasp_recommended_pp import VASP_RECOMMENDED_PP
from drm_copt_estefania.drm_copt_estefania_tools import main_path, incar_tags_bulk


materials = {
    'Co': 'mp-102',
    'Pt': 'mp-126',
    #  'CoPt': 'mp-949'   # problem: it is the primitive cell, prepare the conventional manually
}

for material in materials:

    atoms = read(f'{main_path}/bulk/{material}.poscar')

    incar_tags = deepcopy(incar_tags_bulk)

    if 'Co' in material:
        incar_tags['ISPIN'] = 2

    job = StructureOptimization(
        atoms=atoms,
        incar_tags=incar_tags,
        kspacing=0.66,
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