"""从命令行验证候选路线检查器的几何口径与退出状态。由 Codex 辅助编写。"""

import json
from pathlib import Path
import subprocess
import sys
import unittest


SCRIPT = Path(__file__).with_name("route_audit.py")


def run_audit(points, *args):
    return subprocess.run(
        [sys.executable, "-X", "utf8", "-B", str(SCRIPT), "-", *args],
        input=json.dumps(points),
        capture_output=True,
        encoding="utf-8",
        check=False,
    )


class RouteAuditTests(unittest.TestCase):
    def get_report(self, points, status=0, *args):
        result = run_audit(points, *args)
        self.assertEqual(result.returncode, status, result.stderr)
        return json.loads(result.stdout)

    def test_direct_route_outside_rounded_corner_is_legal(self):
        report = self.get_report([[20.8, 14.8], [20, 14]])
        self.assertAlmostEqual(report["length_m"], 1.1313708499)
        self.assertAlmostEqual(report["min_building_gap_m"], 1.6970562748)
        self.assertEqual(report["violations"], [])

    def test_legal_endpoints_do_not_make_the_corner_cut_legal(self):
        report = self.get_report([[20.5, 16], [22, 14.5]], 1)
        self.assertAlmostEqual(report["min_building_gap_m"], 1.0606601718)
        self.assertEqual(report["violations"], ["building_clearance"])
        self.assertEqual(report["segments"][0]["nearest_building"], "B1")

    def test_route_may_touch_but_not_enter_building(self):
        tangent = self.get_report([[20, 14.5], [37, 14.5]])
        self.assertAlmostEqual(tangent["min_building_gap_m"], 1.5)
        crossing = self.get_report([[20, 20], [37, 20]], 1)
        self.assertEqual(crossing["min_building_gap_m"], 0)

    def test_budget_and_field_boundary_are_inclusive(self):
        path = [[1.5, 2], [1.5, 92]]
        report = self.get_report(path, 0, "--max-seconds", "30")
        self.assertEqual(report["travel_seconds"], 30)
        self.assertEqual(report["min_field_gap_m"], 1.5)
        over = self.get_report(path, 1, "--max-seconds", "29.99")
        self.assertEqual(over["violations"], ["travel_time"])
        outside = self.get_report([[1.49, 2], [1.5, 92]], 1)
        self.assertEqual(outside["violations"], ["field_clearance"])

    def test_requested_endpoint_cannot_be_skipped(self):
        report = self.get_report([[20, 14]], 1, "--end", "21", "14")
        self.assertEqual(report["violations"], ["end_mismatch"])
        self.assertEqual(report["travel_seconds"], 0)

    def test_polyline_checks_later_segments(self):
        report = self.get_report([[20, 12], [20, 15], [36, 20]], 1)
        self.assertGreaterEqual(report["segments"][0]["building_gap_m"], 1.5)
        self.assertLess(report["segments"][1]["building_gap_m"], 1.5)
        self.assertEqual(report["violations"], ["building_clearance"])

    def test_nonfinite_coordinate_is_an_input_error(self):
        result = run_audit([[0, float("nan")]])
        self.assertEqual(result.returncode, 2)
        self.assertIn("坐标", result.stderr)


if __name__ == "__main__":
    unittest.main()
