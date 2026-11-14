# Babymodel: monomer → oligomer → plaque execution spec

## Biology (poster language)
- **Start:** neurons, inactivated_microglia, activated_microglia only (no initial plaques or senescent cells).
- Neurons secrete **amyloid_monomer** (APP/Aβ proxy).
- **Monomer → oligomer** via voxel kinetics in `custom.cpp` when local monomer exceeds a threshold.
- **Oligomer → amyloid_plaque** cell agent via nucleation near neurons; oligomer/monomer mass is debited from the field.
- Microglia **uptake** monomer (fastest) and oligomer (slower); **phagocytose** plaques (slowest).
- **Senescent_microglia** arise only by rules from **activated** microglia under chronic oligomer/cytokine (immobile, limited clearance).
- Plaques do **not** secrete Aβ; they can **grow** by consuming local oligomer (reaccumulation after partial clearance).

## PhysiCell files touched
- `config/PhysiCell_settings.xml` — substrates, neuron cell type, secretion/uptake/chemotaxis
- `config/cell_rules.csv` — CBHG v3, **8 columns**, no header
- `config/cells.csv` — add `neuron` positions
- `custom_modules/custom.cpp` — conversion, nucleation, plaque growth, colors, zero plaque secretion

## Build & run
```bash
cd /Users/genevievelab/Documents/PhysiCell
make load PROJ=babymodel
make -j6
rm -rf output && mkdir output
./project config/PhysiCell_settings.xml
```

## Tunables (`user_parameters` + custom.cpp)
- Monomer→oligomer threshold and rate
- Plaque nucleation oligomer threshold and probability
- Neuron secretion rate / target
