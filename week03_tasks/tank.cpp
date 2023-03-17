// Task 3

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <functional>

double bisection_one_dim(
    const std::function<double(double)> &f,
    double left,
    double right,
    double epsilon_f,
    double epsilon_a
) {
    while (fabs(right - left) >= epsilon_a && fabs(f(left)) >= epsilon_f) {
        auto mid = (right + left) / 2;
        auto val = f(mid);
        if (val <= 0) {
            left = mid;
        }
        else {
            right = mid;
        }
    }

    return left;
}

namespace tank {
    struct Point3D {
        double x, y, z;
    };

    Point3D operator- (const Point3D& p1, const Point3D& p2) {
        return Point3D{p1.x - p2.x, p1.y - p2.y, p1.z - p2.z};
    }

    double dot(const Point3D& p1, const Point3D& p2) {
        return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
    }

    Point3D vec(const Point3D& p1, const Point3D& p2) {
        return Point3D{p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x};
    }

    struct Triangle {
        Point3D a, b, c;
    };

    double vol(const Triangle& t, const Point3D& d) {
        return dot(t.a - d, vec(t.b - d, t.c - d)) / 6.0;
    }

    class TankSkeleton {
    public:
        void load(const std::string& path) {
            std::ifstream fin(path);
            std::string s;
            std::vector<Point3D> pts;
            bool first_point = true;
            while (getline(fin, s)) {
                if (s.find("endsolid OpenSCAD_Model") != std::string::npos) {
                    break;
                }
                if (s.find("vertex") == std::string::npos && !pts.empty()) {
                    triangles.push_back(Triangle{pts[0], pts[1], pts[2]});
                    pts.clear();
                }
                else if (s.find("vertex") != std::string::npos) {
                    std::istringstream ss(s);
                    std::string vertex;
                    double x, y, z;
                    ss >> vertex >> x >> y >> z;
                    if (first_point) {
                        min_z = z;
                        max_z = z;
                        first_point = false;
                    }
                    else {
                        min_z = std::min(min_z, z);
                        max_z = std::max(max_z, z);
                    }
                    pts.push_back(Point3D{x, y, z});
                }
            }
            fin.close();
        }

        double water_volume(double z) {
            Point3D d{0, 0, z};
            double volume = 0;
            for (auto& t: triangles) {
                if (t.a.z > z && t.b.z > z && t.c.z > z) {
                    continue;
                }
                volume += vol(t, d);
            }

            return volume;
        }

        double water_level(double volume) {
            std::function<double(double)> f = [this, volume](double z)
                    { return water_volume(z) - volume; };
            double left = min_z, right = max_z;
            double eps_f = 0.00001, eps_a = 0.000001;
            return bisection_one_dim(f, left, right, eps_f, eps_a);
        }
    private:
        std::vector<Triangle> triangles;
        double min_z = 0.0, max_z = 0.0;
    };
}

int main() {
    tank::TankSkeleton ts;
    ts.load("/home/dword/Desktop/Files/reps/NumericMethods/week03_tasks/tank.stl");
    std::cout << ts.water_volume(ts.water_level(5000)) << std::endl;
    return 0;
}

