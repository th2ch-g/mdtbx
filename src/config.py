from pymol import cmd

# pymol setting
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

# tleap settings
SYSTEM_NAME = "SYS"
AVOGADRO_CONST = 6.022  # 10**23 mol^-1

# ref: https://pubs.acs.org/doi/10.1021/acs.jcim.1c00794
WATER_DENSITY = 0.997  # g/cm^3
WATER_WEGITH = 18  # g/mol
WATER_VOLUME = WATER_WEGITH / WATER_DENSITY / AVOGADRO_CONST / 100  # nm^3/number

TIP3P_DENSITY = 0.980  # g/cm^3
TIP3P_VOLUME = WATER_WEGITH / TIP3P_DENSITY / AVOGADRO_CONST / 100  # nm^3/number

TIP4P_DENSITY = 0.994  # g/cm^3
TIP4P_VOLUME = WATER_WEGITH / TIP4P_DENSITY / AVOGADRO_CONST / 100  # nm^3/number

TIP5P_DENSITY = 0.985  # g/cm^3
TIP5P_VOLUME = WATER_WEGITH / TIP5P_DENSITY / AVOGADRO_CONST / 100  # nm^3/number

OPC_DENSITY = 0.997  # g/cm^3
OPC_VOLUME = WATER_WEGITH / OPC_DENSITY / AVOGADRO_CONST / 100  # nm^3/number

# gaussian setting
# ref: https://qiita.com/Ag_smith/items/430e9efb32a855d4c511
GAUSSIAN_CMD = "g16"
STRUCTURE_OPTIMIZATION = "#p opt=(tight) scf=qc b3lyp/6-31+g(d,p)"
SINGLE_POINT_CALCULATION = "#p hf/6-31g(d) pop=mk iop(6/33=2,6/42=6) scf=tight"
