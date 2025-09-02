from pymol import cmd


def select_protein(model=None):
    common = "polymer.protein or (resn ACE) or (resn NME)"
    if model is None:
        cmd.select("prot", common)
    else:
        cmd.select("prot", f"{common} & {model}")


cmd.extend("select_protein", select_protein)


def select_water(model=None):
    common = "(resn TIP3P or resn HOH or resn OPC or resn TIP3 or resn WAT or resn TP3)"
    if model is None:
        cmd.select("water", common)
    else:
        cmd.select("water", f"{common} & {model}")


cmd.extend("select_water", select_water)


def select_ion(model=None):
    common = "(resn K\\+ or resn Na\\+ or resn Cl- or resn CL or resn Cl or resn K or resn Na or resn NA)"
    if model is None:
        cmd.select("ion", common)
    else:
        cmd.select("ion", f"{common} & {model}")


cmd.extend("select_ion", select_ion)


def select_lipid(model=None):
    common = "(resn PC or resn PA or resn OL or resn CHL)"
    if model is None:
        cmd.select("lipid", common)
    else:
        cmd.select("lipid", f"{common} & {model}")


cmd.extend("select_lipid", select_lipid)


def select_nonh(model=None):
    common = "not name H*"
    if model is None:
        cmd.select("nonh", common)
    else:
        cmd.select("nonh", f"{common} & {model}")


cmd.extend("select_nonh", select_nonh)


def select_template():
    select_protein()
    select_nonh()
    select_lipid()
    select_water()
    select_ion()
    cmd.hide("everything", "ion or water or lipid")
    cmd.select("None")


cmd.extend("select_template", select_template)


def find_ssbond(target_resname: str = "CYS", target_atomname: str = "SG", cutoff=3.0):
    selection = f"resn {target_resname} and name {target_atomname}"
    pairs = cmd.find_pairs(selection, selection, cutoff=cutoff)
    already_set = set()
    ssbonds = []
    for idx1, idx2 in pairs:
        if idx1[1] in already_set or idx2[1] in already_set:
            continue
        already_set.add(idx1[1])
        already_set.add(idx2[1])
        tmp = []
        cmd.iterate_state(1, f"index {idx1[1]}", "tmp.append(resi)", space=locals())
        cmd.iterate_state(1, f"index {idx2[1]}", "tmp.append(resi)", space=locals())
        ssbonds.append(tmp)
    ssbond_str_lines = []
    for res1, res2 in ssbonds:
        ssbond_str_lines.append(
            f"bond SYS.{res1}.{target_atomname} SYS.{res2}.{target_atomname}"
        )
    ssbond_str = "\n".join(ssbond_str_lines)
    print("amber ssbond commands:")
    print(ssbond_str)


cmd.extend("find_ssbond", find_ssbond)
