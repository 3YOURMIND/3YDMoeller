# 3YDMoeller #

A C++ single header implementation of Moeller's Triangle-Triangle & Triangle-Box intersection algorithms for BVH trasversal

### Basic Usage ###
Define your vector3 class: it must be an array-like structure implementing the subscript operator and containing 3 floating points, which represent the x,y,z coordinates value. Any of the following would work:
```
std::vector<float>
std::array<double,3> 
long double[3]
```

The provided methods are:

```c++
threeyd::moeller::TriangleIntersects<myVecType>::box(myVecType& Va, myVecType& Vb, myVecType& Vc, myVecType& boxCenter, myVecType& boxHalfSize);
```

```c++
threeyd::moeller::TriangleIntersects<myVecType>::triangle(myVecType& Va, myVecType& Vb, myVecType& Vc, myVecType& Ua, myVecType& Ub, myVecType& Uc)
```

If you want the Endpoints of the intersection line between two triangles, just pass them to the triangle-triangle intersection test. The boolean reference keeps track of coplanarity.

```c++
threeyd::moeller::TriangleIntersects<myVecType>::triangle(myVecType& Va, myVecType& Vb, myVecType& Vc, myVecType& Ua, myVecType& Ub, myVecType& Uc, myVecType& outIntersectionLineEndPoint1, myVecType& outIntersectionLineEndPoint2, bool& outIsCoplanar)
```

### References ###

Oller, Tomas. (2004). A Fast Triangle-Triangle Intersection Test. Journal of Graphic Tools. 2. 10.1080/10867651.1997.10487472. 