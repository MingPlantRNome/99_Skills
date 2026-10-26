"""Visual authoring/software contract only, never a biological validity test."""
import copy
import importlib.util
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("report_renderer", ROOT / "scripts/render_report_figures.py")
renderer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(renderer)
DATA = json.loads((ROOT / "templates/report-figures.example.json").read_text())


def test_manifest_and_svg_generation(tmp_path):
    index = renderer.generate(DATA, tmp_path)
    assert len(index["figures"]) == 7
    assert index["purpose"] == "template"
    for f in index["figures"]:
        svg = (tmp_path / f["files"][0]).read_text()
        assert ET.fromstring(svg).tag.endswith("svg")
        assert "格式示例｜非真实设计" in svg
        assert "非比例、坐标未定" in svg
    assert json.loads((tmp_path / "figure-index.json").read_text()) == index


def test_hypothesis_line_and_evidence_text():
    svg = renderer.render_figure(DATA["figures"][0], "template")
    assert 'stroke-dasharray="7 5"' in svg and "待验证" in svg
    assert DATA["figures"][0]["links"][0]["basis"]


@pytest.mark.parametrize("fid", ["../escape", "/absolute", "..", "has/slash", "x.svg"])
def test_filename_cannot_escape_output(fid, tmp_path):
    data = copy.deepcopy(DATA)
    data["figures"][0]["id"] = fid
    with pytest.raises(ValueError, match="figure id"):
        renderer.generate(data, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_duplicate_ids_rejected():
    data = copy.deepcopy(DATA)
    data["figures"][1]["id"] = data["figures"][0]["id"]
    with pytest.raises(ValueError, match="duplicate"):
        renderer.validate(data)


@pytest.mark.parametrize("status", ["direct", "indirect", "computed"])
def test_evidence_requires_identifiable_source(status):
    data = copy.deepcopy(DATA)
    data["figures"][0]["links"][0]["status"] = status
    with pytest.raises(ValueError, match="source"):
        renderer.validate(data)
    data["figures"][0]["links"][0]["source"] = "local-assay-01"
    renderer.validate(data)


def test_unsupported_status_and_dense_graph_rejected():
    data = copy.deepcopy(DATA)
    data["figures"][0]["links"][0]["status"] = "guaranteed"
    with pytest.raises(ValueError):
        renderer.validate(data)
    data = copy.deepcopy(DATA)
    data["figures"][0]["nodes"] *= 2
    with pytest.raises(ValueError, match="split dense"):
        renderer.validate(data)


def test_label_is_xml_escaped():
    data = copy.deepcopy(DATA)
    data["figures"][0]["nodes"][0]["label"] = '<script>" &'
    renderer.validate(data)
    svg = renderer.render_figure(data["figures"][0], "project")
    assert "<script>" not in svg and "&lt;script&gt;" in svg and "&amp;" in svg
    ET.fromstring(svg)


def test_packaged_template_images_all_exist():
    report = (ROOT / "templates/design-report.md").read_text()
    images = re.findall(r"!\[[^]]*\]\(([^)]+)\)", report)
    assert len(images) == 7
    for image in images:
        assert (ROOT / "templates" / image).is_file()
        assert (ROOT / "templates" / image).with_suffix(".svg").is_file()
    assert "非比例、坐标未定" in report


def test_hub_routes_visual_style_and_preserves_boundary_contract():
    hub = (ROOT / "SKILL.md").read_text()
    rules = (ROOT / "references/report-visual-style.md").read_text()
    assert "references/report-visual-style.md" in hub
    assert "templates/report-figures.example.json" in hub
    assert "scripts/render_report_figures.py" in hub
    for term in ["不能只交付代码块", "替代解释", "坐标未定", "原始 URL", "不能表示证据强度"]:
        assert term in rules


def test_newline_rejected_to_bound_label_layout():
    data = copy.deepcopy(DATA)
    data["figures"][0]["nodes"][0]["label"] = "label\nsecond line"
    with pytest.raises(ValueError, match="control"):
        renderer.validate(data)
