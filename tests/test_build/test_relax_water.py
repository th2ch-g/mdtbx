import argparse
import types

import numpy as np
import pytest

from src.build import relax_water


def _parse_args(argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    relax_water.add_subcmd(subparsers)
    return parser.parse_args(argv)


def test_defaults():
    args = _parse_args(
        [
            "relax_water",
            "-p",
            "leap.parm7",
            "-x",
            "leap.rst7",
        ]
    )

    assert args.output == "relaxed.rst7"
    assert args.pdb_output == "relaxed.pdb"
    assert args.steps == 2000
    assert args.minimize_iterations == 100
    assert args.water_resname == ["WAT", "HOH", "SOL"]


def test_box_lengths_and_angles_for_orthogonal_box():
    box = relax_water._box_lengths_and_angles(
        np.array(
            [
                [30.0, 0.0, 0.0],
                [0.0, 40.0, 0.0],
                [0.0, 0.0, 50.0],
            ]
        )
    )

    assert box == pytest.approx([30.0, 40.0, 50.0, 90.0, 90.0, 90.0])


def test_add_solute_restraints_excludes_water():
    from openmm import System, Vec3, unit
    from openmm.app import Topology, element

    topology = Topology()
    chain = topology.addChain()
    solute = topology.addResidue("ALA", chain)
    water = topology.addResidue("WAT", chain)
    topology.addAtom("CA", element.carbon, solute)
    topology.addAtom("O", element.oxygen, water)
    positions = [Vec3(0, 0, 0), Vec3(1, 1, 1)] * unit.nanometer
    system = System()
    system.addParticle(12.0)
    system.addParticle(16.0)

    restrained = relax_water._add_solute_restraints(
        system,
        topology,
        positions,
        ["WAT"],
        10.0,
    )

    assert restrained == 1
    assert system.getNumForces() == 1
    assert system.getForce(0).getNumParticles() == 1


def test_run_validates_and_calls_relax(tmp_path, monkeypatch):
    parm = tmp_path / "leap.parm7"
    rst = tmp_path / "leap.rst7"
    parm.touch()
    rst.touch()
    calls = []
    monkeypatch.setattr(relax_water, "_relax", lambda args: calls.append(args))
    args = types.SimpleNamespace(
        parm=str(parm),
        rst=str(rst),
        output=str(tmp_path / "relaxed.rst7"),
        pdb_output=str(tmp_path / "relaxed.pdb"),
        steps=2000,
        minimize_iterations=100,
        temperature=300.0,
        restraint=10.0,
        water_resname=["WAT"],
        seed=2026,
    )

    relax_water.run(args)

    assert calls == [args]


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("steps", 0, "--steps must be positive"),
        (
            "minimize_iterations",
            -1,
            "--minimize-iterations must be non-negative",
        ),
        ("temperature", 0, "--temperature must be positive"),
        ("restraint", -1, "--restraint must be non-negative"),
    ],
)
def test_invalid_numeric_argument(tmp_path, field, value, message):
    parm = tmp_path / "leap.parm7"
    rst = tmp_path / "leap.rst7"
    parm.touch()
    rst.touch()
    args = types.SimpleNamespace(
        parm=str(parm),
        rst=str(rst),
        steps=2000,
        minimize_iterations=100,
        temperature=300.0,
        restraint=10.0,
    )
    setattr(args, field, value)

    with pytest.raises(ValueError, match=message):
        relax_water._validate_args(args)
