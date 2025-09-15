from pymol import cmd


def f(pdb_id: str):
    cmd.fetch(pdb_id)


cmd.extend("f", f)
