from pymol import cmd
from pathlib import Path


def color_plddt(model=None):
    if model is None:
        cmd.color("0x0053D6", "b < 100")
        cmd.color("0x65CBF3", "b < 90")
        cmd.color("0xFFDB13", "b < 70")
        cmd.color("0xFF7D45", "b < 50")
    else:
        cmd.color("0x0053D6", f"{model} and b < 100")
        cmd.color("0x65CBF3", f"{model} and b < 90")
        cmd.color("0xFFDB13", f"{model} and b < 70")
        cmd.color("0xFF7D45", f"{model} and b < 50")


cmd.extend("color_plddt", color_plddt)


def color_baker():
    cmd.spectrum("count", "lightorange_lightpink_lightblue", "polymer.protein")


cmd.extend("color_baker", color_baker)


def ray_png(png="out_pymol_ray.png"):
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.ray(500, 500)
    cmd.png(png, dpi=300)
    print(f"ray image is saved as {png}")
    cmd.bg_color("black")


cmd.extend("ray_png", ray_png)


def load_as(pdb_files: str, model_name="target"):
    print(f"loading {pdb_files} as {model_name}")
    for path in Path(".").glob(pdb_files):
        print(path)
        cmd.load(path, model_name)


cmd.extend("load_as", load_as)
