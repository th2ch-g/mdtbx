from pymol import cmd
from .selector import *  # NoQA
from .builder import *  # NoQA
from .visualizer import *  # NoQA
from .arrow import *  # NoQA


# similar to alignto
# but super has a different algorithm from alignto
def super_all():
    obj_list = cmd.get_object_list("all")
    for obj in obj_list[1:]:
        cmd.super(obj, obj_list[0])


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


def mydefault():
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


cmd.extend("mydefault", mydefault)
mydefault()

print("Loadded!")
