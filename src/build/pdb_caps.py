def normalize_methyl_hydrogen_names(line: str) -> str:
    """Normalize PyMOL methyl hydrogen names to AMBER-style PDB names."""
    for index in range(1, 4):
        replacement = f" H{index} "
        line = line.replace(f"{index}HH3", replacement, 1)
        line = line.replace(f"HH3{index}", replacement, 1)
    return line
