from copy import deepcopy
from ase.build import surface
from ase.io import read
from ase.visualize import view
from ase.build.tools import sort
from ase.constraints import FixAtoms
from vasptools.structureopt import StructureOptimization
from vasptools.vasp_recommended_pp import VASP_RECOMMENDED_PP
from drm_copt_estefania.drm_copt_estefania_tools import main_path, incar_tags_slab


for material in ['Co', 'Pt', 'CoPt']:

    bulk = read(f"{main_path}/bulk/{material}/CONTCAR")

    # Create the (111) surface and expand to a 2x2 supercell
    slab = surface(bulk, (1, 1, 1), layers=4, vacuum=8, periodic=True)
    slab = slab.repeat((2, 2, 1))

    # Move the atoms away from the boundaries of the supercell
    slab.translate((0.5, 0.5, 0.0))
    slab.pbc = (True, True, False)
    slab.wrap()

    # Constrain the bottom two layers
    fix_constraint = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < 0.5*slab.cell[2, 2]])
    slab.set_constraint(fix_constraint)

    incar_tags = deepcopy(incar_tags_slab)
    if 'Co' in material:
        incar_tags['ISPIN'] = 2
        for key in ['LDIPOL', 'IDIPOL', 'DIPOL']:
            incar_tags_slab.pop(key, None)

    slab = sort(slab)
    job = StructureOptimization(
        atoms=slab,
        incar_tags=incar_tags,
        kspacing=1.0,
        kpointstype='gamma',
        potcar_dict=VASP_RECOMMENDED_PP,
        periodicity='2d',
        magmom={'Co': 1.67}
    )
    view(slab)
    #job.write_input_files(folder_name=f"{main_path}/slabs/{material}_111")
