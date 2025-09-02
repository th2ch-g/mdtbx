from collections import defaultdict
from pymol import cmd, editor


def addace(selection: str = None):
    if selection is None:
        selection = ""
    else:
        selection = " and " + selection
    editor.attach_amino_acid(f"first (polymer.protein {selection}) and name N", "ace")
    print("ACE added")
    cmd.rebuild()


cmd.extend("addace", addace)

"""
when saving the pdb file, run the following command
sed -i -e 's/1HH3 ACE/ H1  ACE/g' b.pdb
sed -i -e 's/2HH3 ACE/ H2  ACE/g' b.pdb
sed -i -e 's/3HH3 ACE/ H3  ACE/g' b.pdb
"""


def addnme(selection: str = None):
    if selection is None:
        selection = ""
    else:
        selection = " and " + selection
    cmd.select("oxts", f"name OXT {selection}")
    cmd.remove("oxts")
    print("OXT atoms removed")
    editor.attach_amino_acid(f"last (polymer.protein {selection}) and name C", "nme")
    print("NME added")


cmd.extend("addnme", addnme)

"""
when saving the pdb file, run the following command
sed -i -e '/3HH3 NME /a TER' b.pdb
sed -i -e 's/CH3 NME/C   NME/g' b.pdb
sed -i -e 's/1HH3 NME/ H1  NME/g' b.pdb
sed -i -e 's/2HH3 NME/ H2  NME/g' b.pdb
sed -i -e 's/3HH3 NME/ H3  NME/g' b.pdb
"""


def rename_atomname_renumber(selection: str):
    element_dict = defaultdict(int)
    for atom in cmd.get_model(selection).atom:
        element = atom.name[0]
        element_dict[element] += 1
        cmd.alter(
            f"{selection} and id {atom.id}", f"name='{element}{element_dict[element]}'"
        )
        print(f"{atom.name} renamed to {element}{element_dict[element]}")
    cmd.rebuild()
    print(f"Final counts: {dict(element_dict)}")


cmd.extend("rename_atomname_renumber", rename_atomname_renumber)
