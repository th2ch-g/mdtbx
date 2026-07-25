import importlib.util
from pathlib import Path


PLUGIN_PATH = (
    Path(__file__).resolve().parents[2]
    / "pymol-plugins"
    / "pymol_plugins"
    / "selector.py"
)


class FakeCmd:
    def __init__(self):
        self.select_calls = []
        self.deselect_calls = 0

    def select(self, name, expression):
        self.select_calls.append((name, expression))

    def extend(self, *_args):
        return None

    def hide(self, *_args):
        return None

    def deselect(self):
        self.deselect_calls += 1


def _load_selector():
    spec = importlib.util.spec_from_file_location("selector_under_test", PLUGIN_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_model_filter_applies_to_entire_protein_expression():
    selector = _load_selector()
    fake_cmd = FakeCmd()
    selector.cmd = fake_cmd

    selector.select_protein("model_a")

    assert fake_cmd.select_calls == [
        (
            "prot",
            "(polymer.protein or (resn ACE) or (resn NME)) and (model_a)",
        )
    ]


def test_select_template_clears_selection_with_deselect():
    selector = _load_selector()
    fake_cmd = FakeCmd()
    selector.cmd = fake_cmd

    selector.select_template()

    assert fake_cmd.deselect_calls == 1
