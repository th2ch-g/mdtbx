from pymol import cmd
from .selector import *  # NoQA
from .builder import *  # NoQA
from .visualizer import *  # NoQA
from .arrow import *  # NoQA
from .alias import *  # NoQA


# similar to alignto
# but super has a different algorithm from alignto
def super_all():
    obj_list = cmd.get_object_list("all")
    for obj in obj_list[1:]:
        cmd.super(obj, obj_list[0])
    cmd.zoom(obj_list[0])


cmd.extend("super_all", super_all)


def ls_res(selection):
    residues = set()
    cmd.iterate_state(1, selection, "residues.add((resi, resn))", space=locals())
    residue_list = list(residues)
    residue_list.sort(key=lambda x: int(x[0]))
    for resi, resn in residue_list:
        print(f" Residue: {resi}  Name: {resn}")


cmd.extend("ls_res", ls_res)


# set translate command
# this is effective when you use laptop without mouse
# select -> extract object -> MovO
def set_movo():
    cmd.button("L", "CtSh", "MovO")


cmd.extend("set_movo", set_movo)


def comdist(selection1, selection2):
    cmd.pseudoatom("com1", selection1)
    cmd.pseudoatom("com2", selection2)
    cmd.distance("comdist1-2", "com1", "com2")


cmd.extend("comdist", comdist)


# use gui for instead
# def make_trj_movie():
#     cmd.viewport(12000, 10000)
#     cmd.mset("1 -%d" % cmd.count_states())
#     cmd.mpng("traj_frame_")
#     import subprocess
#
#     subprocess.run(
#         "ffmpeg -r 30 -i traj_frame_%04d.png -c:v libx264 -pix_fmt yuv420p out.mov",
#         shell=True,
#         check=True,
#     )
#     subprocess.run("rm -f traj_frame_*.png", shell=True, check=True)
#     print("movie is saved as out.mp4")
#
#
# cmd.extend("make_trj_movie", make_trj_movie)


def default_settings():
    cmd.reinitialize("settings")


cmd.extend("default_settings", default_settings)


def laboo_settings():
    default_settings()
    cmd.set("sphere_scale", 0.22)
    cmd.set("sphere_scale", 0.22, "elem C+N+O+S+Cl+F+Na+Mg")
    cmd.set("sphere_scale", 0.13, "elem H")
    cmd.set("hide_long_bonds", 1)
    cmd.set("dash_gap", 0)
    cmd.set("dash_gap", 0.15)
    cmd.set("dash_length", 0.05)
    cmd.set("dash_round_ends", 0)
    cmd.set("dash_radius", 0.05)
    cmd.set("label_size", 18)
    cmd.set("cartoon_loop_radius", 0.1)
    cmd.set("cartoon_putty_radius", 0.2)
    cmd.set("cartoon_oval_length", 0.8)
    cmd.set("cartoon_gap_cutoff", 0)
    cmd.set("cartoon_rect_length", 1.2)
    cmd.set("label_digits", 3)
    cmd.set("ray_opaque_background", 0)
    cmd.set("dash_length", 0.2500)
    cmd.set("retain_order", 1)
    cmd.set("pdb_retain_ids", 1)
    cmd.set("suspend_undo_atom_count", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.05)
    cmd.set("ray_trace_color", "gray37")
    # cmd.bg_color("white")
    cmd.set("fog", 1)
    cmd.set("ambient", 0.66)
    cmd.set("reflect", 0)
    cmd.set("spec_reflect", 0.08)
    cmd.set("light_count", 2)


cmd.extend("laboo_settings", laboo_settings)


def my_settings():
    default_settings()
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


cmd.extend("my_settings", my_settings)
my_settings()
