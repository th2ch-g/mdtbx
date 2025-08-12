from pymol import cmd

# save setting
cmd.set("pdb_reformat_names_mode", 2)
cmd.set("retain_order", 1)
cmd.set("pdb_retain_ids", 1)
cmd.set("pdb_conect_all", "off")
# cmd.set("cartoon_highlight_color", "grey50")

# visualization setting
cmd.set("ray_shadows", "0")
cmd.set("ray_trace_mode", "1")
cmd.set("spec_reflect", "0")
cmd.set("cartoon_tube_radius", "0.5")
cmd.set("ambient", "0.4")
cmd.set("ray_trace_color", "black")

# grompp setting
MAXWARN = 10
