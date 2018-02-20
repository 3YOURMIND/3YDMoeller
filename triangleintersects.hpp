#ifndef TRIANGLEINTERSECTS_H
#define TRIANGLEINTERSECTS_H
#include <utility>
#include <type_traits>
#include <iostream>

#include <math.h>
#include <stdio.h>

/**
 * Single header mplementation of Triangle-Triangle and Triangle-Box intersection tests by Tomas Akenine Moeller
 * @see http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt
 * @see http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tritri_isectline.txt
*/

namespace moeller
{
template <typename...>
using void_t = void;

template <class T, class Index, typename = void>
struct has_subscript_operator : std::false_type
{
};

template <class T, class Index>
struct has_subscript_operator<T, Index, void_t<decltype(std::declval<T>()[std::declval<Index>()])>> : std::true_type
{
};

/**
 * Template @param T is any array like structure of any floating point type that implements the subscript operator
 *
 * Provides bool moeller:TriangleIntersects<T>::triangle(T v1, T v2, T v3, T u1, T u2, T u3);
 * Provides bool moeller:TriangleIntersects<T>::triangle(T v1, T v2, T v3, T u1, T u2, T u3, T out_inters_endpoint1,
 *  T out_inters_endpoint1, bool out_is_coplanar);
 * Provides bool moeller:TriangleIntersects<T>::box(T v1, T v2, T v3, T boxCenter, T boxHalfSize);
 */

template <class TemplatedVec>
class TriangleIntersects
{
    static_assert(has_subscript_operator<TemplatedVec, size_t>::value, "Vector3 needs to implement [] operator");
    typedef typename std::decay<decltype(std::declval<TemplatedVec>()[0])>::type declfloat;
    static_assert(std::is_floating_point<declfloat>::value, "The elements of vector3 need to be floats or doubles");

  public:
    static bool triangle(TemplatedVec &firstV1, TemplatedVec &firstV2, TemplatedVec &firstV3, TemplatedVec &secondV1,
                         TemplatedVec &secondV2, TemplatedVec &secondV3)
    {
        TemplatedVec IntersectionLineEndPoint1;
        TemplatedVec IntersectionLineEndPoint2;
        bool coplanar;
        return tri_tri_intersect_with_isectline(firstV1, firstV2, firstV3, secondV1, secondV2, secondV3, coplanar,
                                                IntersectionLineEndPoint1, IntersectionLineEndPoint2, false);
    }

    static bool triangle(TemplatedVec &firstV1, TemplatedVec &firstV2, TemplatedVec &firstV3, TemplatedVec &secondV1,
                         TemplatedVec &secondV2, TemplatedVec &secondV3, TemplatedVec &IntersectionLineEndPoint1,
                         TemplatedVec &IntersectionLineEndPoint2, bool &coplanar)
    {
        return tri_tri_intersect_with_isectline(firstV1, firstV2, firstV3, secondV1, secondV2, secondV3, coplanar,
                                                IntersectionLineEndPoint1, IntersectionLineEndPoint2, true);
    }

    static bool box(TemplatedVec triangleV1, TemplatedVec triangleV2, TemplatedVec triangleV3, TemplatedVec boxCenter,
                    TemplatedVec boxHalfSize)
    {
        return triBoxOverlap(triangleV1, triangleV2, triangleV3, boxCenter, boxHalfSize);
    }

  private:
    // Constants definitions
    static constexpr declfloat EPSILON = 0.000001f;
    static const size_t X = 0;
    static const size_t Y = 1;
    static const size_t Z = 2;

    // Helper methods
    inline static void cross(TemplatedVec &dest, const TemplatedVec &v1, const TemplatedVec &v2)
    {
        dest[0] = v1[1] * v2[2] - v1[2] * v2[1];
        dest[1] = v1[2] * v2[0] - v1[0] * v2[2];
        dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
    }
    inline static declfloat dot(const TemplatedVec &v1, const TemplatedVec &v2)
    {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }
    inline static void sub(TemplatedVec &dest, const TemplatedVec &v1, const TemplatedVec &v2)
    {
        dest[0] = v1[0] - v2[0];
        dest[1] = v1[1] - v2[1];
        dest[2] = v1[2] - v2[2];
    }
    inline static void add(TemplatedVec &dest, const TemplatedVec &v1, const TemplatedVec &v2)
    {
        dest[0] = v1[0] + v2[0];
        dest[1] = v1[1] + v2[1];
        dest[2] = v1[2] + v2[2];
    }
    inline static void mult(TemplatedVec &dest, const TemplatedVec &v, const declfloat factor)
    {
        dest[0] = factor * v[0];
        dest[1] = factor * v[1];
        dest[2] = factor * v[2];
    }
    inline static void set(TemplatedVec &dest, const TemplatedVec &src)
    {
        dest[0] = src[0];
        dest[1] = src[1];
        dest[2] = src[2];
    }
    inline static void sort(declfloat &a, declfloat &b)
    {
        if (a > b)
        {
            std::swap(a, b);
        }
    }
    inline static bool sort2(declfloat &a, declfloat &b)
    {
        if (a > b)
        {
            std::swap(a, b);
            return true;
        }
        return false;
    }

    inline static void findMinMax(const declfloat x0, const declfloat x1, const declfloat x2, declfloat &min,
                                  declfloat &max)
    {
        min = std::min(std::min(x0, x1), x2);
        max = std::max(std::max(x0, x1), x2);
    }

    // Tests for Box-Triangle
    inline static bool axisTestX01(const TemplatedVec &v0, const TemplatedVec &v2, const TemplatedVec &boxhalfsize,
                                   const declfloat a, const declfloat b, const declfloat fa, const declfloat fb,
                                   declfloat &min, declfloat &max, declfloat &rad)
    {
        declfloat p0 = a * v0[Y] - b * v0[Z];
        declfloat p2 = a * v2[Y] - b * v2[Z];
        if (p0 < p2)
        {
            min = p0;
            max = p2;
        }
        else
        {
            min = p2;
            max = p0;
        }
        rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];
        if (min > rad || max < -rad)
        {
            return false;
        }
        return true;
    }
    inline static bool axisTestX2(const TemplatedVec &v0, const TemplatedVec &v1, const TemplatedVec &boxhalfsize,
                                  const declfloat a, const declfloat b, const declfloat fa, const declfloat fb,
                                  declfloat &min, declfloat &max, declfloat &rad)
    {
        declfloat p0 = a * v0[Y] - b * v0[Z];
        declfloat p1 = a * v1[Y] - b * v1[Z];
        if (p0 < p1)
        {
            min = p0;
            max = p1;
        }
        else
        {
            min = p1;
            max = p0;
        }
        rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];
        if (min > rad || max < -rad)
        {
            return false;
        }
        return true;
    }
    inline static bool axisTestY02(const TemplatedVec &v0, const TemplatedVec &v2, const TemplatedVec &boxhalfsize,
                                   const declfloat a, const declfloat b, const declfloat fa, const declfloat fb,
                                   declfloat &min, declfloat &max, declfloat &rad)
    {
        declfloat p0 = -a * v0[X] + b * v0[Z];
        declfloat p2 = -a * v2[X] + b * v2[Z];
        if (p0 < p2)
        {
            min = p0;
            max = p2;
        }
        else
        {
            min = p2;
            max = p0;
        }
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];
        if (min > rad || max < -rad)
        {
            return false;
        }
        return true;
    }
    inline static bool axisTestY1(const TemplatedVec &v0, const TemplatedVec &v1, const TemplatedVec &boxhalfsize,
                                  const declfloat a, const declfloat b, const declfloat fa, const declfloat fb,
                                  declfloat &min, declfloat &max, declfloat &rad)
    {
        declfloat p0 = -a * v0[X] + b * v0[Z];
        declfloat p1 = -a * v1[X] + b * v1[Z];
        if (p0 < p1)
        {
            min = p0;
            max = p1;
        }
        else
        {
            min = p1;
            max = p0;
        }
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];
        if (min > rad || max < -rad)
        {
            return false;
        }
        return true;
    }
    inline static bool axisTestZ12(const TemplatedVec &v1, const TemplatedVec &v2, const TemplatedVec &boxhalfsize,
                                   const declfloat a, const declfloat b, const declfloat fa, const declfloat fb,
                                   declfloat &min, declfloat &max, declfloat &rad)
    {
        declfloat p1 = a * v1[X] - b * v1[Y];
        declfloat p2 = a * v2[X] - b * v2[Y];
        if (p2 < p1)
        {
            min = p2;
            max = p1;
        }
        else
        {
            min = p1;
            max = p2;
        }
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];
        if (min > rad || max < -rad)
        {
            return false;
        }
        return true;
    }
    inline static bool axisTestZ0(const TemplatedVec &v0, const TemplatedVec &v1, const TemplatedVec &boxhalfsize,
                                  const declfloat a, const declfloat b, const declfloat fa, const declfloat fb,
                                  declfloat &min, declfloat &max, declfloat &rad)
    {
        declfloat p0 = a * v0[X] - b * v0[Y];
        declfloat p1 = a * v1[X] - b * v1[Y];
        if (p0 < p1)
        {
            min = p0;
            max = p1;
        }
        else
        {
            min = p1;
            max = p0;
        }
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];
        if (min > rad || max < -rad)
        {
            return false;
        }
        return true;
    }

    // Tests for Triangle-Triangle
    inline static bool edge_edge_test(const TemplatedVec &V0, const TemplatedVec &U0, const TemplatedVec &U1,
                                      const size_t i0, const size_t i1, declfloat Ax, declfloat Ay)
    {
        declfloat Bx, By, Cx, Cy, f, d, e;
        Bx = U0[i0] - U1[i0];
        By = U0[i1] - U1[i1];
        Cx = V0[i0] - U0[i0];
        Cy = V0[i1] - U0[i1];
        f = Ay * Bx - Ax * By;
        d = By * Cx - Bx * Cy;
        if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f))
        {
            e = Ax * Cy - Ay * Cx;
            if (f > 0)
            {
                if (e >= 0 && e <= f) return true;
            }
            else
            {
                if (e <= 0 && e >= f) return true;
            }
        }
        return false;
    }
    inline static bool edge_against_tri_edge(const TemplatedVec &V0, const TemplatedVec &V1, const TemplatedVec &U0,
                                             const TemplatedVec &U1, const TemplatedVec &U2, const size_t i0,
                                             const size_t i1)
    {
        declfloat Ax, Ay;
        Ax = V1[i0] - V0[i0];
        Ay = V1[i1] - V0[i1];
        /* test edge U0,U1 against V0,V1 */
        if (edge_edge_test(V0, U0, U1, i0, i1, Ax, Ay))
        {
            return true;
        }
        /* test edge U1,U2 against V0,V1 */
        if (edge_edge_test(V0, U1, U2, i0, i1, Ax, Ay))
        {
            return true;
        }
        /* test edge U2,U1 against V0,V1 */
        if (edge_edge_test(V0, U2, U0, i0, i1, Ax, Ay))
        {
            return true;
        }
        return false;
    }

    inline static void isect2(const TemplatedVec &VTX0, const TemplatedVec &VTX1, const TemplatedVec &VTX2,
                              const declfloat VV0, const declfloat VV1, const declfloat VV2, const declfloat D0,
                              const declfloat D1, const declfloat D2, declfloat &isect0, declfloat &isect1,
                              TemplatedVec &isectpoint0, TemplatedVec &isectpoint1)
    {
        declfloat tmp = D0 / (D0 - D1);
        TemplatedVec diff;
        isect0 = VV0 + (VV1 - VV0) * tmp;
        sub(diff, VTX1, VTX0);
        mult(diff, diff, tmp);
        add(isectpoint0, diff, VTX0);
        tmp = D0 / (D0 - D2);
        isect1 = VV0 + (VV2 - VV0) * tmp;
        sub(diff, VTX2, VTX0);
        mult(diff, diff, tmp);
        add(isectpoint1, VTX0, diff);
    }

    inline static bool point_in_tri(const TemplatedVec &V0, const TemplatedVec &U0, const TemplatedVec &U1,
                                    const TemplatedVec &U2, const size_t i0, const size_t i1)
    {
        declfloat a, b, c, d0, d1, d2;
        /* is T1 completly inside T2? */
        /* check if V0 is inside tri(U0,U1,U2) */
        a = U1[i1] - U0[i1];
        b = -(U1[i0] - U0[i0]);
        c = -a * U0[i0] - b * U0[i1];
        d0 = a * V0[i0] + b * V0[i1] + c;

        a = U2[i1] - U1[i1];
        b = -(U2[i0] - U1[i0]);
        c = -a * U1[i0] - b * U1[i1];
        d1 = a * V0[i0] + b * V0[i1] + c;

        a = U0[i1] - U2[i1];
        b = -(U0[i0] - U2[i0]);
        c = -a * U2[i0] - b * U2[i1];
        d2 = a * V0[i0] + b * V0[i1] + c;
        if (d0 * d1 > 0.0)
        {
            if (d0 * d2 > 0.0) return true;
        }
        else
            return false;
    }

    // Private methods
    static bool planeBoxOverlap(const TemplatedVec &normal, const TemplatedVec &vert,
                                const TemplatedVec &maxbox)  // -NJMP-
    {
        size_t q;
        declfloat v;
        TemplatedVec vmin, vmax;
        for (q = X; q <= Z; q++)
        {
            v = vert[q];
            if (normal[q] > 0.0f)
            {
                vmin[q] = -maxbox[q] - v;
                vmax[q] = maxbox[q] - v;
            }
            else
            {
                vmin[q] = maxbox[q] - v;
                vmax[q] = -maxbox[q] - v;
            }
        }
        if (dot(normal, vmin) > 0.0f)
        {
            return false;
        }
        if (dot(normal, vmax) >= 0.0f)
        {
            return true;
        }
        return false;
    }

    static bool triBoxOverlap(const TemplatedVec &trivert0, const TemplatedVec &trivert1, const TemplatedVec &trivert2,
                              const TemplatedVec &boxcenter, const TemplatedVec &boxhalfsize)
    {
        /*    use separating axis theorem to test overlap between triangle and box */
        /*    need to test for overlap in these directions: */
        /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
        /*       we do not even need to test these) */
        /*    2) normal of the triangle */
        /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
        /*       this gives 3x3=9 more tests */

        TemplatedVec v0, v1, v2;
        declfloat min, max, rad, fex, fey, fez;
        TemplatedVec normal, e0, e1, e2;
        /* This is the fastest branch on Sun */
        /* move everything so that the boxcenter is in (0,0,0) */
        sub(v0, trivert0, boxcenter);
        sub(v1, trivert1, boxcenter);
        sub(v2, trivert2, boxcenter);
        /* compute triangle edges */
        sub(e0, v1, v0); /* tri edge 0 */
        sub(e1, v2, v1); /* tri edge 1 */
        sub(e2, v0, v2); /* tri edge 2 */

        /* Bullet 3:  */
        /*  test the 9 tests first (this was faster) */
        fex = fabsf(e0[X]);
        fey = fabsf(e0[Y]);
        fez = fabsf(e0[Z]);

        if (!axisTestX01(v0, v2, boxhalfsize, e0[Z], e0[Y], fez, fey, min, max, rad))
        {
            return false;
        }
        if (!axisTestY02(v0, v2, boxhalfsize, e0[Z], e0[X], fez, fex, min, max, rad))
        {
            return false;
        }
        if (!axisTestZ12(v1, v2, boxhalfsize, e0[Y], e0[X], fey, fex, min, max, rad))
        {
            return false;
        }

        fex = fabsf(e1[X]);
        fey = fabsf(e1[Y]);
        fez = fabsf(e1[Z]);

        if (!axisTestX01(v0, v2, boxhalfsize, e1[Z], e1[Y], fez, fey, min, max, rad))
        {
            return false;
        }
        if (!axisTestY02(v0, v2, boxhalfsize, e1[Z], e1[X], fez, fex, min, max, rad))
        {
            return false;
        }
        if (!axisTestZ0(v0, v1, boxhalfsize, e1[Y], e1[X], fey, fex, min, max, rad))
        {
            return false;
        }

        fex = fabsf(e2[X]);
        fey = fabsf(e2[Y]);
        fez = fabsf(e2[Z]);

        if (!axisTestX2(v0, v1, boxhalfsize, e2[Z], e2[Y], fez, fey, min, max, rad))
        {
            return false;
        }
        if (!axisTestY1(v0, v1, boxhalfsize, e2[Z], e2[X], fez, fex, min, max, rad))
        {
            return false;
        }
        if (!axisTestZ12(v1, v2, boxhalfsize, e2[Y], e2[X], fey, fex, min, max, rad))
        {
            return false;
        }

        /* Bullet 1: */
        /*  first test overlap in the {x,y,z}-directions */
        /*  find min, max of the triangle each direction, and test for overlap in */
        /*  that direction -- this is equivalent to testing a minimal AABB around */
        /*  the triangle against the AABB */
        /* test in X-direction */

        findMinMax(v0[X], v1[X], v2[X], min, max);
        if (min > boxhalfsize[X] || max < -boxhalfsize[X])
        {
            return false;
        }
        /* test in Y-direction */
        findMinMax(v0[Y], v1[Y], v2[Y], min, max);
        if (min > boxhalfsize[Y] || max < -boxhalfsize[Y])
        {
            return false;
        }
        /* test in Z-direction */
        findMinMax(v0[Z], v1[Z], v2[Z], min, max);
        if (min > boxhalfsize[Z] || max < -boxhalfsize[Z])
        {
            return false;
        }

        /* Bullet 2: */
        /*  test if the box intersects the plane of the triangle */
        /*  compute plane equation of triangle: normal*x+d=0 */
        cross(normal, e0, e1);

        if (!planeBoxOverlap(normal, v0, boxhalfsize))
        {
            return false;
        }
        return true; /* box and triangle overlaps */
    }

    static bool tri_tri_intersect_with_isectline(const TemplatedVec &V0, const TemplatedVec &V1, const TemplatedVec &V2,
                                                 const TemplatedVec &U0, const TemplatedVec &U1, const TemplatedVec &U2,
                                                 bool &coplanar, TemplatedVec &isectpt1, TemplatedVec &isectpt2,
                                                 bool check_isect_endpoints)
    {
        TemplatedVec E1 = {0.f, 0.f, 0.f};
        TemplatedVec E2 = {0.f, 0.f, 0.f};
        TemplatedVec N1 = {0.f, 0.f, 0.f};
        TemplatedVec N2 = {0.f, 0.f, 0.f};
        declfloat d1, d2;
        declfloat du0, du1, du2, dv0, dv1, dv2;
        TemplatedVec D = {0.f, 0.f, 0.f};
        declfloat isect1[2], isect2[2];
        TemplatedVec isectpointA1 = {0.f, 0.f, 0.f};
        TemplatedVec isectpointA2 = {0.f, 0.f, 0.f};
        TemplatedVec isectpointB1 = {0.f, 0.f, 0.f};
        TemplatedVec isectpointB2 = {0.f, 0.f, 0.f};
        declfloat du0du1, du0du2, dv0dv1, dv0dv2;
        size_t index;
        declfloat vp0, vp1, vp2;
        declfloat up0, up1, up2;
        declfloat b, c, max;
        size_t smallest1, smallest2;

        /* compute plane equation of triangle(V0,V1,V2) */
        sub(E1, V1, V0);
        sub(E2, V2, V0);
        cross(N1, E1, E2);
        d1 = -dot(N1, V0);
        /* plane equation 1: N1.X+d1=0 */

        /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
        du0 = dot(N1, U0) + d1;
        du1 = dot(N1, U1) + d1;
        du2 = dot(N1, U2) + d1;

        if (fabs(du0) < EPSILON) du0 = 0.0;
        if (fabs(du1) < EPSILON) du1 = 0.0;
        if (fabs(du2) < EPSILON) du2 = 0.0;

        du0du1 = du0 * du1;
        du0du2 = du0 * du2;

        if (du0du1 > 0.0f && du0du2 > 0.0f) /* same sign on all of them + not equal 0 ? */
            return false;                   /* no intersection occurs */

        /* compute plane of triangle (U0,U1,U2) */
        sub(E1, U1, U0);
        sub(E2, U2, U0);
        cross(N2, E1, E2);
        d2 = -dot(N2, U0);
        /* plane equation 2: N2.X+d2=0 */

        /* put V0,V1,V2 into plane equation 2 */
        dv0 = dot(N2, V0) + d2;
        dv1 = dot(N2, V1) + d2;
        dv2 = dot(N2, V2) + d2;

        if (fabs(dv0) < EPSILON) dv0 = 0.0;
        if (fabs(dv1) < EPSILON) dv1 = 0.0;
        if (fabs(dv2) < EPSILON) dv2 = 0.0;

        dv0dv1 = dv0 * dv1;
        dv0dv2 = dv0 * dv2;

        if (dv0dv1 > 0.0f && dv0dv2 > 0.0f) /* same sign on all of them + not equal 0 ? */
            return false;                   /* no intersection occurs */

        /* compute direction of intersection line */
        cross(D, N1, N2);

        /* compute and index to the largest component of D */
        max = fabs(D[0]);
        index = 0;
        b = fabs(D[1]);
        c = fabs(D[2]);
        if (b > max) max = b, index = 1;
        if (c > max) max = c, index = 2;

        /* this is the simplified projection onto L*/
        vp0 = V0[index];
        vp1 = V1[index];
        vp2 = V2[index];

        up0 = U0[index];
        up1 = U1[index];
        up2 = U2[index];

        /* compute interval for triangle 1 */
        coplanar = compute_intervals_isectline(V0, V1, V2, vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, isect1[0],
                                               isect1[1], isectpointA1, isectpointA2);
        if (coplanar)
        {
            return coplanar_tri_tri(N1, V0, V1, V2, U0, U1, U2);
        }

        /* compute interval for triangle 2 */
        compute_intervals_isectline(U0, U1, U2, up0, up1, up2, du0, du1, du2, du0du1, du0du2, isect2[0], isect2[1],
                                    isectpointB1, isectpointB2);

        smallest1 = sort2(isect1[0], isect1[1]);
        smallest2 = sort2(isect2[0], isect2[1]);

        if (isect1[1] < isect2[0] || isect2[1] < isect1[0]) return false;

        /* at this point, we know that the triangles intersect */

        if (!check_isect_endpoints)
        {
            return true;
        }

        if (isect2[0] < isect1[0])
        {
            if (smallest1 == 0)
            {
                dot(isectpt1, isectpointA1);
            }
            else
            {
                dot(isectpt1, isectpointA2);
            }

            if (isect2[1] < isect1[1])
            {
                if (smallest2 == 0)
                {
                    dot(isectpt2, isectpointB2);
                }
                else
                {
                    dot(isectpt2, isectpointB1);
                }
            }
            else
            {
                if (smallest1 == 0)
                {
                    dot(isectpt2, isectpointA2);
                }
                else
                {
                    dot(isectpt2, isectpointA1);
                }
            }
        }
        else
        {
            if (smallest2 == 0)
            {
                dot(isectpt1, isectpointB1);
            }
            else
            {
                dot(isectpt1, isectpointB2);
            }

            if (isect2[1] > isect1[1])
            {
                if (smallest1 == 0)
                {
                    dot(isectpt2, isectpointA2);
                }
                else
                {
                    dot(isectpt2, isectpointA1);
                }
            }
            else
            {
                if (smallest2 == 0)
                {
                    dot(isectpt2, isectpointB2);
                }
                else
                {
                    dot(isectpt2, isectpointB1);
                }
            }
        }
        return true;
    }

    static bool coplanar_tri_tri(const TemplatedVec &N, const TemplatedVec &V0, const TemplatedVec &V1,
                                 const TemplatedVec &V2, const TemplatedVec &U0, const TemplatedVec &U1,
                                 const TemplatedVec &U2)
    {
        TemplatedVec A = {0.f, 0.f, 0.f};
        size_t i0, i1;
        /* first project onto an axis-aligned plane, that maximizes the area */
        /* of the triangles, compute indices: i0,i1. */
        A[0] = fabs(N[0]);
        A[1] = fabs(N[1]);
        A[2] = fabs(N[2]);
        if (A[0] > A[1])
        {
            if (A[0] > A[2])
            {
                i0 = 1; /* A[0] is greatest */
                i1 = 2;
            }
            else
            {
                i0 = 0; /* A[2] is greatest */
                i1 = 1;
            }
        }
        else /* A[0]<=A[1] */
        {
            if (A[2] > A[1])
            {
                i0 = 0; /* A[2] is greatest */
                i1 = 1;
            }
            else
            {
                i0 = 0; /* A[1] is greatest */
                i1 = 2;
            }
        }

        /* test all edges of triangle 1 against the edges of triangle 2 */
        if (edge_against_tri_edge(V0, V1, U0, U1, U2, i0, i1))
        {
            return true;
        }
        if (edge_against_tri_edge(V1, V2, U0, U1, U2, i0, i1))
        {
            return true;
        }
        if (edge_against_tri_edge(V2, V0, U0, U1, U2, i0, i1))
        {
            return true;
        }

        /* finally, test if tri1 is totally contained in tri2 or vice versa */
        if (point_in_tri(V0, U0, U1, U2, i0, i1))
        {
            return true;
        }
        if (point_in_tri(U0, V0, V1, V2, i0, i1))
        {
            return true;
        }

        return false;
    }

    inline static bool compute_intervals_isectline(const TemplatedVec &VERT0, const TemplatedVec &VERT1,
                                                   const TemplatedVec &VERT2, const declfloat VV0, const declfloat VV1,
                                                   const declfloat VV2, const declfloat D0, const declfloat D1,
                                                   const declfloat D2, const declfloat D0D1, const declfloat D0D2,
                                                   declfloat &isect0, declfloat &isect1, TemplatedVec &isectpoint0,
                                                   TemplatedVec &isectpoint1)
    {
        if (D0D1 > 0.0f)
        {
            /* here we know that D0D2<=0.0 */
            /* that is D0, D1 are on the same side, D2 on the other or on the plane */
            isect2(VERT2, VERT0, VERT1, VV2, VV0, VV1, D2, D0, D1, isect0, isect1, isectpoint0, isectpoint1);
        }
        else if (D0D2 > 0.0f)
        {
            /* here we know that d0d1<=0.0 */
            isect2(VERT1, VERT0, VERT2, VV1, VV0, VV2, D1, D0, D2, isect0, isect1, isectpoint0, isectpoint1);
        }
        else if (D1 * D2 > 0.0f || D0 != 0.0f)
        {
            /* here we know that d0d1<=0.0 or that D0!=0.0 */
            isect2(VERT0, VERT1, VERT2, VV0, VV1, VV2, D0, D1, D2, isect0, isect1, isectpoint0, isectpoint1);
        }
        else if (D1 != 0.0f)
        {
            isect2(VERT1, VERT0, VERT2, VV1, VV0, VV2, D1, D0, D2, isect0, isect1, isectpoint0, isectpoint1);
        }
        else if (D2 != 0.0f)
        {
            isect2(VERT2, VERT0, VERT1, VV2, VV0, VV1, D2, D0, D1, isect0, isect1, isectpoint0, isectpoint1);
        }
        else
        {
            /* triangles are coplanar */
            return true;
        }
        return false;
    }
};
}
#endif  // TRIANGLEINTERSECTS_H
