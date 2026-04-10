'''
PyMOL Morphing Command Plugin

Linear coordinate interpolation for Open-Source PyMOL.
'''

import numpy as np


def morph(sele1, sele2, name=None, state1=1, state2=1,
          steps=30, superpose=0, _self=None):
    '''
DESCRIPTION

    Generate a morphing trajectory between two conformations using linear
    coordinate interpolation. Works with Open-Source PyMOL.

USAGE

    morph sele1, sele2 [, name [, state1 [, state2 [, steps [, superpose ]]]]]

ARGUMENTS

    sele1     = str: start conformation object/selection
    sele2     = str: end conformation object/selection
    name      = str: output object name {default: auto}
    state1    = int: state of sele1 {default: 1}
    state2    = int: state of sele2 {default: 1}
    steps     = int: number of output states {default: 30}
    superpose = 0/1: align sele1 onto sele2 before morphing {default: 0}

EXAMPLE

    morph 1ake, 4ake
    morph open, closed, name=transition, steps=50, superpose=1
    '''
    if _self is None:
        from pymol import cmd as _self

    if not name:
        name = _self.get_unused_name('morph')

    state1 = int(state1)
    state2 = int(state2)
    steps = int(steps)
    superpose = int(superpose)

    if superpose:
        _self.align(sele1, sele2, mobile_state=state1, target_state=state2)

    coords1 = _self.get_coords(sele1, state1)
    coords2 = _self.get_coords(sele2, state2)

    if coords1 is None or coords2 is None:
        print(' morph-error: failed to get coordinates')
        return

    if coords1.shape != coords2.shape:
        print(f' morph-error: atom count mismatch ({coords1.shape[0]} vs {coords2.shape[0]})')
        return

    # Build multi-state object: all states share topology from sele1
    for i in range(1, steps + 1):
        _self.create(name, sele1, state1, i)

    # Overwrite each state with linearly interpolated coordinates
    for i in range(steps):
        t = i / (steps - 1) if steps > 1 else 0.0
        coords_i = (1.0 - t) * coords1 + t * coords2
        _self.load_coords(coords_i, name, state=i + 1)

    print(f' morph: created "{name}" with {steps} states')


def __init_plugin__(app=None):
    from pymol import cmd
    cmd.extend('morph', morph)
