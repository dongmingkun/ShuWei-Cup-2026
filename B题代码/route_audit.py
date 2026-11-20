"""核验给定的圆盘车辆折线路线。由 Codex（OpenAI）辅助编写，使用日期 2026-11-20。"""

import argparse
import json
import math
from pathlib import Path
import sys

from geometry import (
    BUILDINGS,
    EPS,
    FIELD_SIZE,
    SPEED,
    VEHICLE_RADIUS,
    dist_point_to_rect,
    segment_hits_box,
)


def _point(raw):
    if not isinstance(raw, (list, tuple)) or len(raw) != 2:
        raise ValueError("每个路点须为 [x, y] 坐标")
    if any(type(x) not in (float, int) for x in raw):
        raise ValueError("坐标须为数值")
    values = tuple(map(float, raw))
    if not all(map(math.isfinite, values)):
        raise ValueError("坐标须为有限数值")
    return values


def _point_to_segment(point, a, b, length):
    if length == 0:
        return math.dist(point, a)
    ux, uy = (b[0] - a[0]) / length, (b[1] - a[1]) / length
    projection = (point[0] - a[0]) * ux + (point[1] - a[1]) * uy
    if projection <= 0:
        return math.dist(point, a)
    if projection >= length:
        return math.dist(point, b)
    return math.hypot(point[0] - (a[0] + projection * ux),
                      point[1] - (a[1] + projection * uy))


def _building_gap(a, b, length, building):
    rect = building.x0, building.x1, building.y0, building.y1
    if segment_hits_box((*a, 0.0), (*b, 0.0), (*rect, 0.0, 0.0), eps=0.0):
        return 0.0
    corners = ((x, y) for x in rect[:2] for y in rect[2:])
    return min(
        dist_point_to_rect(*a, *rect),
        dist_point_to_rect(*b, *rect),
        *(_point_to_segment(corner, a, b, length) for corner in corners),
    )


def _field_gap(point):
    x, y = point
    return min(x, y, FIELD_SIZE - x, FIELD_SIZE - y)


def audit_route(points, *, max_seconds=None, start=None, end=None):
    """逐段核对净距、端点与时间预算；返回报告，不寻找或修改路线。"""
    if not isinstance(points, (list, tuple)) or not points:
        raise ValueError("路线须为非空的坐标数组")
    points = [_point(p) for p in points]
    if max_seconds is not None and (not math.isfinite(max_seconds) or max_seconds < 0):
        raise ValueError("时间上限须为有限非负数")
    segments = []
    field_illegal = building_illegal = False
    legs = zip(points, points[1:]) if len(points) > 1 else [(points[0], points[0])]
    for index, (a, b) in enumerate(legs, start=1):
        length = math.dist(a, b)
        if not math.isfinite(length):
            raise ValueError("路段长度超出浮点表示范围")
        gap, building_index = min(
            (_building_gap(a, b, length, obstacle), i)
            for i, obstacle in enumerate(BUILDINGS, start=1)
        )
        field_gap = min(_field_gap(a), _field_gap(b))
        field_illegal |= field_gap < VEHICLE_RADIUS - EPS
        building_illegal |= gap < VEHICLE_RADIUS - EPS
        segments.append({
            "index": index,
            "length_m": length,
            "field_gap_m": field_gap,
            "building_gap_m": gap,
            "nearest_building": f"B{building_index}",
        })
    length = math.fsum(segment["length_m"] for segment in segments)
    travel_seconds = length / SPEED
    violations = []
    if field_illegal:
        violations.append("field_clearance")
    if building_illegal:
        violations.append("building_clearance")
    if start is not None and math.dist(_point(start), points[0]) > EPS:
        violations.append("start_mismatch")
    if end is not None and math.dist(_point(end), points[-1]) > EPS:
        violations.append("end_mismatch")
    if max_seconds is not None and travel_seconds > max_seconds + EPS:
        violations.append("travel_time")
    return {
        "valid_route": not violations,
        "length_m": length,
        "travel_seconds": travel_seconds,
        "min_field_gap_m": min(segment["field_gap_m"] for segment in segments),
        "min_building_gap_m": min(segment["building_gap_m"] for segment in segments),
        "violations": violations,
        "segments": segments,
    }


def main():
    parser = argparse.ArgumentParser(description="验算给定折线路线的圆盘净距和行驶时间")
    parser.add_argument("route", help="UTF-8 JSON 文件，或 - 从标准输入读取 [[x,y], ...]")
    parser.add_argument("--max-seconds", type=float, help="行驶时间上限（秒）")
    parser.add_argument("--start", nargs=2, type=float, metavar=("X", "Y"), help="预期起点")
    parser.add_argument("--end", nargs=2, type=float, metavar=("X", "Y"), help="预期终点")
    args = parser.parse_args()
    try:
        raw = sys.stdin.read() if args.route == "-" else Path(args.route).read_text(encoding="utf-8-sig")
        result = audit_route(
            json.loads(raw), max_seconds=args.max_seconds, start=args.start, end=args.end,
        )
        report = json.dumps(result, indent=2, allow_nan=False)
    except (OSError, ValueError, OverflowError) as error:
        parser.error(str(error))
    print(report)
    return 0 if result["valid_route"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
