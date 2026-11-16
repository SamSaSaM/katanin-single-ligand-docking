#!/usr/bin/env python3
"""Analyze heavy-atom contacts and plot the best Vina pose in 3D.

Distance contacts are descriptive only and do not establish hydrogen bonds or
experimental binding-site membership.
"""

import csv
from math import dist
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
NAME = "R1_2ET_R2_4OE4design_R3_amide"
RECEPTOR = ROOT / "data/receptor/katanin.pdbqt"
POSE = ROOT / "results" / f"{NAME}_docked.pdbqt"
RESULTS = ROOT / "results"
CUTOFF = 4.0

ELEMENT_COLORS = {
    "C": "#4D4D4D",
    "N": "#2457FF",
    "O": "#E53935",
    "S": "#E6B800",
    "P": "#F28E2B",
}
COVALENT_RADII = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07}


def read_atoms(path, first_model=False):
    records = []
    for line in Path(path).read_text().splitlines():
        if first_model and line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        name = line[12:16].strip()
        residue = line[17:20].strip()
        chain = line[21:22].strip()
        number = line[22:26].strip()
        xyz = tuple(float(line[index:index + 8]) for index in (30, 38, 46))
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = "".join(char for char in name if char.isalpha())[:1].upper()
        # Normalize common AutoDock atom types to chemical elements.
        element = {
            "A": "C",
            "OA": "O",
            "NA": "N",
            "SA": "S",
            "HD": "H",
        }.get(element, element)
        if element in {"H", "D"} or name.upper().startswith("H"):
            continue
        records.append(
            {
                "atom": name,
                "residue": residue,
                "chain": chain,
                "number": number,
                "element": element,
                "xyz": xyz,
            }
        )
    return records


def inferred_ligand_bonds(atoms):
    bonds = []
    for i, atom_i in enumerate(atoms):
        for j in range(i + 1, len(atoms)):
            atom_j = atoms[j]
            radius_i = COVALENT_RADII.get(atom_i["element"], 0.77)
            radius_j = COVALENT_RADII.get(atom_j["element"], 0.77)
            distance = dist(atom_i["xyz"], atom_j["xyz"])
            if 0.4 < distance <= 1.25 * (radius_i + radius_j):
                bonds.append((i, j))
    return bonds


def main():
    receptor_atoms = read_atoms(RECEPTOR)
    ligand_atoms = read_atoms(POSE, first_model=True)
    if not receptor_atoms or not ligand_atoms:
        raise SystemExit("Could not parse receptor or ligand coordinates")

    closest_by_residue = {}
    for receptor_atom in receptor_atoms:
        residue_key = (
            receptor_atom["residue"],
            receptor_atom["chain"],
            receptor_atom["number"],
        )
        for ligand_atom in ligand_atoms:
            separation = dist(receptor_atom["xyz"], ligand_atom["xyz"])
            previous = closest_by_residue.get(residue_key)
            if previous is None or separation < previous["distance_A"]:
                closest_by_residue[residue_key] = {
                    "residue": receptor_atom["residue"],
                    "chain": receptor_atom["chain"],
                    "residue_number": receptor_atom["number"],
                    "distance_A": separation,
                    "receptor_atom": receptor_atom["atom"],
                    "receptor_element": receptor_atom["element"],
                    "receptor_xyz": receptor_atom["xyz"],
                    "ligand_atom": ligand_atom["atom"],
                    "ligand_element": ligand_atom["element"],
                    "ligand_xyz": ligand_atom["xyz"],
                }

    contacts = sorted(
        (
            contact
            for contact in closest_by_residue.values()
            if contact["distance_A"] <= CUTOFF
        ),
        key=lambda contact: contact["distance_A"],
    )

    csv_path = RESULTS / "contact_summary.csv"
    with csv_path.open("w", newline="") as handle:
        fields = [
            "residue",
            "chain",
            "residue_number",
            "distance_A",
            "receptor_atom",
            "receptor_element",
            "ligand_atom",
            "ligand_element",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for contact in contacts:
            writer.writerow(
                {
                    field: (
                        f"{contact[field]:.2f}"
                        if field == "distance_A"
                        else contact[field]
                    )
                    for field in fields
                }
            )

    figure = plt.figure(figsize=(12, 10), facecolor="white")
    axis = figure.add_subplot(111, projection="3d")

    for first, second in inferred_ligand_bonds(ligand_atoms):
        coordinates = np.array(
            [ligand_atoms[first]["xyz"], ligand_atoms[second]["xyz"]]
        )
        axis.plot(
            coordinates[:, 0],
            coordinates[:, 1],
            coordinates[:, 2],
            color="#222222",
            linewidth=2.5,
            zorder=5,
        )

    for element in sorted({atom["element"] for atom in ligand_atoms}):
        selected = np.array(
            [atom["xyz"] for atom in ligand_atoms if atom["element"] == element]
        )
        axis.scatter(
            selected[:, 0],
            selected[:, 1],
            selected[:, 2],
            s=95,
            color=ELEMENT_COLORS.get(element, "#999999"),
            edgecolor="white",
            linewidth=0.7,
            label=f"Ligand {element}",
            depthshade=True,
            zorder=6,
        )

    for contact in contacts:
        receptor_xyz = np.array(contact["receptor_xyz"])
        ligand_xyz = np.array(contact["ligand_xyz"])
        axis.scatter(
            *receptor_xyz,
            s=85,
            color="#65A765",
            alpha=0.85,
            edgecolor="#2E6B2E",
            linewidth=0.6,
            zorder=3,
        )
        axis.plot(
            [receptor_xyz[0], ligand_xyz[0]],
            [receptor_xyz[1], ligand_xyz[1]],
            [receptor_xyz[2], ligand_xyz[2]],
            linestyle="--",
            color="#9A9A9A",
            linewidth=0.8,
            alpha=0.7,
            zorder=2,
        )
        label = f"{contact['residue']}{contact['residue_number']}"
        axis.text(
            receptor_xyz[0],
            receptor_xyz[1],
            receptor_xyz[2],
            label,
            fontsize=8,
            color="#174D17",
        )

    ligand_xyz = np.array([atom["xyz"] for atom in ligand_atoms])
    center = ligand_xyz.mean(axis=0)
    span = max(np.ptp(ligand_xyz, axis=0).max(), 10.0) / 2 + 3.0
    axis.set_xlim(center[0] - span, center[0] + span)
    axis.set_ylim(center[1] - span, center[1] + span)
    axis.set_zlim(center[2] - span, center[2] + span)
    axis.set_box_aspect((1, 1, 1))
    axis.view_init(elev=22, azim=38)
    axis.set_xlabel("X (Angstrom)")
    axis.set_ylabel("Y (Angstrom)")
    axis.set_zlabel("Z (Angstrom)")
    axis.set_title(
        f"{NAME}\nBest Vina pose: -7.091 kcal/mol; residues within 4.0 Angstrom",
        fontsize=14,
        pad=20,
    )
    axis.grid(False)
    figure.text(
        0.5,
        0.02,
        "Dashed lines show nearest heavy-atom distances, not validated hydrogen bonds.",
        ha="center",
        fontsize=10,
        color="#555555",
    )
    figure.tight_layout(rect=(0, 0.04, 1, 1))
    output_figure = RESULTS / f"{NAME}_binding_pose_3D.png"
    figure.savefig(output_figure, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(f"Saved {csv_path}")
    print(f"Saved {output_figure}")
    print(f"Pocket residues within {CUTOFF:.1f} A: {len(contacts)}")


if __name__ == "__main__":
    main()

