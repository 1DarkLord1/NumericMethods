#include <vector>
#include <list>
#include <cassert>
#include <iostream>

namespace comp_geom {
    struct Point {
        double x, y;
    };

    struct Triangle {
        Point a, b, c;
    };

    Point operator-(const Point &p1, const Point &p2) {
        return Point{p1.x - p2.x, p1.y - p2.y};
    }

    bool operator==(const Point &p1, const Point &p2) {
        return p1.x == p2.x && p1.y == p2.y;
    }

    bool operator!=(const Point &p1, const Point &p2) {
        return p1.x != p2.x || p1.y != p2.y;
    }

    double vec_product(const Point &v1, const Point &v2) {
        return v1.x * v2.y - v2.x * v1.y;
    }

    std::list<Point>::const_iterator get_next(const std::list<Point>& polygon, std::list<Point>::const_iterator cur) {
        cur++;
        if (cur == polygon.end()) {
            cur = polygon.begin();
        }
        return cur;
    }

    std::list<Point>::const_iterator get_prev(const std::list<Point>& polygon, std::list<Point>::const_iterator cur) {
        if (cur == polygon.begin()) {
            cur = --polygon.end();
        }
        else {
            cur--;
        }
        return cur;
    }

    int point_order(const std::list<Point> &polygon, Point p) {
        int order = 0;
        for (auto cur = polygon.begin(); cur != polygon.end(); cur++) {
            auto next = get_next(polygon, cur);

            double t = vec_product(*cur - p, *next - p);
            if (t > 0 && cur->y <= p.y && p.y < next->y) {
                order++;
            } else if (t < 0 && next->y <= p.y && p.y < cur->y) {
                order--;
            }
        }

        return order;
    }

    bool is_ear(const std::list<Point>& polygon, Point prev, Point cur, Point next) {
        std::list<Point> triangle = {prev, cur, next};
        for (auto p: polygon) {
            if (p == prev || p == cur || p == next) {
                continue;
            }
            int p_order = point_order(triangle, p);
            if (p_order != 0) {
                return false;
            }
        }

        return true;
    }

    bool is_convex(Point prev, Point cur, Point next) {
        return vec_product(cur - prev, next - prev) > 0;
    }

    std::vector<Triangle> triangulation(std::list<Point> polygon) {
        const int n = (int)polygon.size();
        std::vector<Triangle> triangles;
        std::list<Point>::const_iterator cur = polygon.begin();

        while (polygon.size() >= 3) {
            auto prev = get_prev(polygon, cur);
            auto next = get_next(polygon, cur);

            Point prev_pt = *prev, cur_pt = *cur, next_pt = *next;

            if (!is_convex(prev_pt, cur_pt, next_pt)) {
                cur++;
                continue;
            }
            if (is_ear(polygon, prev_pt, cur_pt, next_pt)) {
                triangles.push_back(Triangle{prev_pt, cur_pt, next_pt});
                auto cur_copy = cur;
                cur = prev;
                polygon.erase(cur_copy);
            }
            else {
                cur++;
            }
        }

        assert(polygon.size() == 2);
        assert(triangles.size() == n - 2);

        return triangles;
    }
}

int main() {
    std::list<comp_geom::Point> polygon = {
            comp_geom::Point{3, -0.5},
            comp_geom::Point{1, 0},
            comp_geom::Point{0, 2},
            comp_geom::Point{-1, 0},
            comp_geom::Point{-3, -0.5},
            comp_geom::Point{-1, -2},
            comp_geom::Point{0, -4},
            comp_geom::Point{1, -2}
    };

    auto outer = comp_geom::Point{-1.5, 0.5};
    auto inner = comp_geom::Point{-1, -0.5};

    assert(comp_geom::point_order(polygon, outer) == 0);
    assert(comp_geom::point_order(polygon, inner) != 0);

    auto triangles = comp_geom::triangulation(polygon);
    std::cout << "Triangulation: " << std::endl;
    for (int i = 0; i < triangles.size(); i++) {
        auto& t = triangles[i];
        std::cout << "Triangle #" << i + 1 << ":" << std::endl;
        std::cout << "Point #1: " << t.a.x << " " << t.a.y << std::endl;
        std::cout << "Point #2: " << t.b.x << " " << t.b.y << std::endl;
        std::cout << "Point #3: " << t.c.x << " " << t.c.y << std::endl;
    }
    return 0;
}