import shlex

from dram_viz.make_product import main


def test_main(test_annotation_path, tmp_path):
    main(shlex.split(f"--annotations {test_annotation_path} --output-dir {tmp_path}/output"), standalone_mode=False)
    assert (tmp_path / "output").exists()
    assert (tmp_path / "output" / "product.html").exists()
    assert (tmp_path / "output" / "product.tsv").exists()
