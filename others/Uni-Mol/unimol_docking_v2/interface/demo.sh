python /Users/flash/Desktop/FlashDock/Uni-Mol/unimol_docking_v2/interface/demo.py --mode single --conf-size 10 --cluster \
        --input-protein 1G9V_receptor.pdb \
        --input-ligand 1G9V_ligand.sdf \
        --input-docking-grid docking_grid.json \
        --output-ligand-name ligand_predict \
        --output-ligand-dir predict_sdf \
        --steric-clash-fix \
        --model-dir /Users/flash/Desktop/Uni-Mol-main/unimol_docking_v2/unimol_docking_v2_240517.pt