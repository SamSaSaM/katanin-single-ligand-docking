# Mode 6: three explicit-H donor-acceptor hydrogen bonds.
reinitialize

load visualization/katanin_receptor.pdb, receptor
load visualization/R1_2ET_R2_4OE4design_R3_amide_mode6.pdb, ligand

hide everything
bg_color white
set ray_opaque_background, off
set antialias, 2

show cartoon, receptor
color gray90, receptor
set cartoon_transparency, 0.80, receptor

select contacts_4A, byres (receptor within 4.0 of ligand)
select pocket_5A, byres (receptor within 5.0 of ligand)
select shell_4to5A, pocket_5A and not contacts_4A
show sticks, pocket_5A
set stick_radius, 0.18, pocket_5A

# Conventional atom colors
color palegreen, contacts_4A and elem C
color wheat, shell_4to5A and elem C
color blue, pocket_5A and elem N
color red, pocket_5A and elem O
color yellow, pocket_5A and elem S
color white, pocket_5A and elem H

show sticks, ligand
set stick_radius, 0.27, ligand
color gray40, ligand and elem C
color blue, ligand and elem N
color red, ligand and elem O
color yellow, ligand and elem S
color white, ligand and elem H

label pocket_5A and name CA, "%s%s" % (resn, resi)
set label_color, black
set label_size, 12
set label_font_id, 7
set label_outline_color, white
set label_position, [1.4, 1.4, 1.4], pocket_5A and name CA
label ligand and elem N+O+S, "%s%s" % (elem, ID)

# Highlight and label the protein N/O atoms directly involved in the three
# explicit-H hydrogen-bond candidates.
select interacting_protein_heavy_atoms, (receptor and chain A and resi 257 and name N) or (receptor and chain A and resi 256 and name OG1) or (receptor and chain A and resi 308 and name OD1)
show spheres, interacting_protein_heavy_atoms
set sphere_scale, 0.32, interacting_protein_heavy_atoms
label receptor and chain A and resi 257 and name N, "LEU257:N"
label receptor and chain A and resi 256 and name OG1, "THR256:OG1"
label receptor and chain A and resi 308 and name OD1, "ASP308:OD1"
set label_color, black, interacting_protein_heavy_atoms
set label_size, 14, interacting_protein_heavy_atoms
set label_outline_color, white, interacting_protein_heavy_atoms
set label_position, [1.0, 1.0, 1.0], interacting_protein_heavy_atoms

# Explicit H...acceptor distances. These pairs passed:
# H...A <= 2.5 A, D...A <= 3.5 A, and D-H...A >= 120 degrees.
# Protein donor LEU257 N-H -> ligand O9
distance hb_LEU257_to_O9, (receptor and chain A and resi 257 and name H), (ligand and id 9)
# Ligand N10-H11 -> THR256 OG1
distance hb_H11_to_THR256, (ligand and id 11), (receptor and chain A and resi 256 and name OG1)
# Ligand N14-H15 -> ASP308 OD1
distance hb_H15_to_ASP308, (ligand and id 15), (receptor and chain A and resi 308 and name OD1)

set dash_color, red, hb_LEU257_to_O9
set dash_color, blue, hb_H11_to_THR256
set dash_color, blue, hb_H15_to_ASP308
set dash_width, 4.0
set dash_gap, 0.22
set dash_length, 0.18
set label_distance_digits, 2
set label_size, 13, hb_LEU257_to_O9
set label_size, 13, hb_H11_to_THR256
set label_size, 13, hb_H15_to_ASP308
set label_color, black
set label_outline_color, white

group explicit_hbonds, hb_LEU257_to_O9 hb_H11_to_THR256 hb_H15_to_ASP308
group mode6_pocket, contacts_4A pocket_5A shell_4to5A interacting_protein_heavy_atoms explicit_hbonds

zoom ligand or pocket_5A, 15
orient ligand
turn y, 18
turn x, 8

ray 2800, 2100
png results/mode6_explicit_H_bonds_5A.png, dpi=300
