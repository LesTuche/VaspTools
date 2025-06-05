import os
from ase.io import read
from vasptools.structureopt import StructureOptimization
from vasptools.vasp_recommended_pp import VASP_RECOMMENDED_PP
from mlp_mxenes_tools import main_path, incar_tags_spe

batch_id = '6'
structures = read(f"{main_path}/spe_calculations/structures_batch_{batch_id}.traj", index=":")
base_output_folder = f'{main_path}/spe_calculations/vasp_files_batch_{batch_id}'

incar_tags_spe['ENCUT'] = 700
incar_tags_spe['NCORE'] = 12  # 12 or 16

os.makedirs(base_output_folder, exist_ok=True)

batch_size = None # e.g. 100, None
for i in range(len(structures)):
    atoms = structures[i]

    if batch_size is None:
        structure_folder = os.path.join(base_output_folder, f"structure_{i + 1}")
    else:
        batch_index = i // batch_size
        batch_folder = os.path.join(base_output_folder, f"batch_{batch_index + 1}")
        os.makedirs(batch_folder, exist_ok=True)
        structure_folder = os.path.join(batch_folder, f"structure_{i + 1}")

    job = StructureOptimization(
        atoms=atoms,
        incar_tags=incar_tags_spe,
        kspacing=0.15,
        kspacing_definition='vasp',
        kpointstype='gamma',
        potcar_dict=VASP_RECOMMENDED_PP,
        periodicity='2d'
    )
    job.write_input_files(folder_name=structure_folder)
