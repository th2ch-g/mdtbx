from pymol import cmd


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
    cmd.ray()
    cmd.png(png, dpi=300)
    print(f"ray image is saved as {png}")

cmd.extend("ray_png", ray_png)
