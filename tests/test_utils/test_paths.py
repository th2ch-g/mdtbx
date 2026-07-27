from mdtbx.utils.paths import ensure_output_parent


def test_ensure_output_parent_creates_nested_directory(tmp_path):
    output_path = tmp_path / "nested" / "output.dat"

    prepared_path = ensure_output_parent(output_path)

    assert prepared_path == output_path
    assert output_path.parent.is_dir()
