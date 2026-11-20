"""对题面口径的回归检查；多警报样例只用于检查时序。由 Codex 辅助编写。"""

import math
import unittest

from demo_smoke import plan
from geometry import clearance_ok, segment_hits_box, shortest_path
from simulator import simulate


class PublicModuleAudit(unittest.TestCase):
    def test_disk_can_leave_southwest_corner_of_b1_directly(self):
        start, goal = (20.8, 14.8), (20.0, 14.0)
        self.assertTrue(clearance_ok(start) and clearance_ok(goal))
        result = shortest_path(start, goal)
        self.assertIsNotNone(result, "整段远离 B1 西南角，净距始终大于 1.5 m")
        length, _ = result
        self.assertAlmostEqual(length, math.dist(start, goal), places=9)

    def test_single_point_contact_with_building_edge_blocks_los(self):
        box = (22.0, 34.0, 16.0, 26.0, 0.0, 20.0)
        self.assertTrue(segment_hits_box((20, 18, 2), (24, 14, 2), box))
        self.assertFalse(segment_hits_box((20, 17.99, 2), (24, 13.99, 2), box))

    def test_reused_vehicle_leaves_after_joint_observation_and_reset(self):
        result = simulate(plan)
        v2 = {(event["alert"], event["kind"]): event for event in result.events[2]}
        observation = v2["T3", "obs"]
        reset = v2["T3", "reset"]
        self.assertAlmostEqual(observation["t_start"], result.alert_start["T3"])
        self.assertAlmostEqual(observation["t_end"] - observation["t_start"], 5.0)
        self.assertAlmostEqual(reset["t_start"], observation["t_end"])
        self.assertGreater(reset["t_end"], reset["t_start"])
        self.assertAlmostEqual(v2["T1", "depart"]["t_start"], reset["t_end"])


if __name__ == "__main__":
    unittest.main()
