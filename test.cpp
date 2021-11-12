#include "triangleintersects.hpp"

#include <array>
#include <cstddef>

struct Vector3
{
    float x;
    float y;
    float z;

    float& operator[](size_t i) noexcept { return *(&x + i); }

    const float& operator[](size_t i) const noexcept { return *(&x + i); }
};

using Triangle = std::array<Vector3, 3>;

int main()
{
    Triangle t0{Vector3{0.f, 1.f, 0.f}, Vector3{1.f, 1.f, 0.f}, Vector3{0.f, -1.f, 0.f}};
    Triangle t1{Vector3{0.f, 1.f, 0.f}, Vector3{1.f, 1.f, 0.f}, Vector3{0.f, -1.f, 0.f}};

    if (!threeyd::moeller::TriangleIntersects<Vector3>::triangle(t0[0], t0[1], t0[2], t1[0], t1[1], t1[2]))
    {
        return 1;
    }

    Vector3 outIntersectionLineEndPoint1;
    Vector3 outIntersectionLineEndPoint2;
    bool outIsCoplanar;
    if (!threeyd::moeller::TriangleIntersects<Vector3>::triangle(t0[0], t0[1], t0[2], t1[0], t1[1], t1[2],
                                                                 outIntersectionLineEndPoint1,
                                                                 outIntersectionLineEndPoint2, outIsCoplanar) ||
        !outIsCoplanar)
    {
        return 1;
    }

    Vector3 boxCenter{0.f, 1.f, 0.f};
    Vector3 boxHalfSize{0.5f, 0.5f, 0.f};
    if (!threeyd::moeller::TriangleIntersects<Vector3>::box(t0[0], t0[1], t0[2], boxCenter, boxHalfSize))
    {
        return 1;
    }

    return 0;
}