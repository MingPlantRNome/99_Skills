"""Check report authoring contract only; not scientific validity of a real report."""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REPORT = (ROOT / "templates/design-report.md").read_text(encoding="utf-8")
RULES = (ROOT / "references/report-boundaries.md").read_text(encoding="utf-8")
HUB = (ROOT / "SKILL.md").read_text(encoding="utf-8")


@pytest.mark.parametrize("heading", [
    "适用范围与不能保证的结论",
    "已知缺点与设计权衡",
    "当前无法回答的问题",
    "待优化方向与验证计划",
])
def test_boundary_topics_are_required_in_template_and_rules(heading):
    assert heading in REPORT and heading in RULES


def test_hub_routes_report_creation_to_boundary_rules():
    assert "references/report-boundaries.md" in HUB
    assert "报告必须包含" in HUB


def test_unknowns_have_reasons_decision_effects_and_validation_endpoints():
    for field in ["无法回答原因", "对本次决策的影响", "最小验证", "可判定结果", "是否阻断交付"]:
        assert field in REPORT


def test_optimization_plan_includes_constraints_controls_and_rejection():
    for field in ["优先级", "触发条件", "保护约束", "证据等级", "风险", "关键对照", "否决判据"]:
        assert field in REPORT
    assert "保留原报道复刻版" in REPORT


def test_missing_sequences_remain_linked_to_machine_blockers():
    assert "unresolved_issues" in REPORT and "unresolved_issues" in RULES
    assert "同步写入" in REPORT and "阻断" in RULES


def test_no_forced_optimization_or_fabricated_benefit_requirement():
    assert "可以建议保留原版" in REPORT
    assert "不填提升倍数" in REPORT
    assert "不能自动证明" in RULES


def test_report_tables_have_consistent_column_counts():
    table = []
    for line in REPORT.splitlines() + [""]:
        if line.startswith("|"):
            table.append(line)
        elif table:
            assert len(table) >= 3
            assert len({row.count("|") for row in table}) == 1
            table = []
