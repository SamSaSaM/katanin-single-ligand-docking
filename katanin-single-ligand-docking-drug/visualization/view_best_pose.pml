# Run from the repository root with: pymol visualization/view_best_pose.pml
load visualization/katanin_receptor.pdb, receptor
load visualization/R1_2ET_R2_4OE4design_R3_amide_mode1.pdb, ligand

hide everything
show cartoon, receptor
color gray80, receptor
set cartoon_transparency, 0.25, receptor

show sticks, ligand
util.cbag ligand
set stick_radius, 0.18, ligand

select pocket, byres (receptor within 4.0 of ligand)
show sticks, pocket
color palegreen, pocket
show sticks, ligand
util.cbag ligand

label pocket and name CA, "%s%s" % (resn, resi)
set label_size, 12
set label_color, black

# Distance-based polar contacts; do not interpret automatically as H-bonds.
distance polar_candidates, (ligand and elem N+O+S), (pocket and elem N+O+S), 3.5
set dash_color, magenta, polar_candidates
hide labels, polar_candidates

bg_color white
set ray_opaque_background, off
zoom ligand, 10
orient ligand
