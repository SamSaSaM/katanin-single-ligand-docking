#!/usr/bin/env python3
"""Reproduce single-ligand docking against prepared KATNA1 receptor."""

import csv
from pathlib import Path
import shutil
import subprocess
import sys

from rdkit import Chem
from rdkit.Chem import AllChem


ROOT = Path(__file__).resolve().parents[1]
NAME = "R1_2ET_R2_4OE4design_R3_amide"
SEED = 20261115
EXHAUSTIVENESS = 32
NUM_MODES = 10


def main():
    for executable in ("vina", "mk_prepare_ligand.py"):
        if shutil.which(executable) is None:
            sys.exit(f"Missing executable: {executable}")

    mol_file = ROOT / "data/ligand" / f"{NAME}.mol"
    receptor = ROOT / "data/receptor/katanin.pdbqt"
    box = ROOT / "config/vina_box.txt"
    output_dir = ROOT / "reproduced_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    mol = Chem.MolFromMolFile(str(mol_file))
    if mol is None:
        sys.exit(f"Could not parse ligand: {mol_file}")

    mol = Chem.AddHs(mol)
    embedding = AllChem.ETKDGv3()
    embedding.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, embedding) != 0:
        sys.exit("Conformer generation failed")

    optimization_status = AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    if optimization_status != 0:
        print("WARNING: MMFF optimization did not fully converge", file=sys.stderr)

    sdf_file = output_dir / f"{NAME}.sdf"
    writer = Chem.SDWriter(str(sdf_file))
    writer.write(mol)
    writer.close()

    ligand_pdbqt = output_dir / f"{NAME}.pdbqt"
    subprocess.run(
        [
            "mk_prepare_ligand.py",
            "-i",
            str(sdf_file),
            "-o",
            str(ligand_pdbqt),
        ],
        check=True,
    )

    docked_file = output_dir / f"{NAME}_docked.pdbqt"
    command = [
        "vina",
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand_pdbqt),
        "--config",
        str(box),
        "--exhaustiveness",
        str(EXHAUSTIVENESS),
        "--num_modes",
        str(NUM_MODES),
        "--seed",
        str(SEED),
        "--out",
        str(docked_file),
    ]
    process = subprocess.run(command, capture_output=True, text=True)
    log_file = output_dir / f"{NAME}_vina.log"
    log_file.write_text(process.stdout + "\n" + process.stderr)
    if process.returncode != 0:
        sys.exit(f"Docking failed; see {log_file}")

    scores = [
        float(line.split()[3])
        for line in docked_file.read_text().splitlines()
        if line.startswith("REMARK VINA RESULT:")
    ]
    summary_file = output_dir / "docking_summary.csv"
    with summary_file.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("compound", "best_vina_kcal_mol", "poses", "output"),
        )
        writer.writeheader()
        writer.writerow(
            {
                "compound": NAME,
                "best_vina_kcal_mol": scores[0] if scores else "",
                "poses": len(scores),
                "output": docked_file.name,
            }
        )

    print(f"Best Vina score: {scores[0] if scores else 'not found'} kcal/mol")
    print(f"Saved results to: {output_dir}")


if __name__ == "__main__":
    main()

