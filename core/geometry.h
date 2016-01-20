
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_GEOMETRY_H
#define PBRT_CORE_GEOMETRY_H

// core/geometry.h*
#include "pbrt.h"

// Geometry Declarations
class Vector {
public:
    // Vector Public Methods
    Vector() { x = y = z = 0.f; }
    Vector(float xx, float yy, float zz)
        : x(xx), y(yy), z(zz) {
        Assert(!HasNaNs());
    }
    bool HasNaNs() const { return isnan(x) || isnan(y) || isnan(z); }
    explicit Vector(const Point &p);
#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    Vector(const Vector &v) {
        Assert(!v.HasNaNs());
        x = v.x; y = v.y; z = v.z;
    }

    Vector &operator=(const Vector &v) {
        Assert(!v.HasNaNs());
        x = v.x; y = v.y; z = v.z;
        return *this;
    }
#endif // !NDEBUG
    Vector operator+(const Vector &v) const {
        Assert(!v.HasNaNs());
        return Vector(x + v.x, y + v.y, z + v.z);
    }

    Vector& operator+=(const Vector &v) {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    Vector operator-(const Vector &v) const {
        Assert(!v.HasNaNs());
        return Vector(x - v.x, y - v.y, z - v.z);
    }

    Vector& operator-=(const Vector &v) {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    Vector operator*(float f) const { return Vector(f*x, f*y, f*z); }

    Vector &operator*=(float f) {
        Assert(!isnan(f));
        x *= f; y *= f; z *= f;
        return *this;
    }
    Vector operator/(float f) const {
        Assert(f != 0);
        float inv = 1.f / f;
        return Vector(x * inv, y * inv, z * inv);
    }

    Vector &operator/=(float f) {
        Assert(f != 0);
        float inv = 1.f / f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    Vector operator-() const { return Vector(-x, -y, -z); }
    float operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    float &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    float LengthSquared() const { return x*x + y*y + z*z; }
    float Length() const { return sqrtf(LengthSquared()); }
    explicit Vector(const Normal &n);

    bool operator==(const Vector &v) const {
        return x == v.x && y == v.y && z == v.z;
    }
    bool operator!=(const Vector &v) const {
        return x != v.x || y != v.y || z != v.z;
    }

    // Vector Public Data
    float x, y, z;
};


class Point {
public:
    // Point Public Methods
    Point() { x = y = z = 0.f; }
    Point(float xx, float yy, float zz)
        : x(xx), y(yy), z(zz) {
        Assert(!HasNaNs());
    }
#ifndef NDEBUG
    Point(const Point &p) {
        Assert(!p.HasNaNs());
        x = p.x; y = p.y; z = p.z;
    }

    Point &operator=(const Point &p) {
        Assert(!p.HasNaNs());
        x = p.x; y = p.y; z = p.z;
        return *this;
    }
#endif // !NDEBUG
    Point operator+(const Vector &v) const {
        Assert(!v.HasNaNs());
        return Point(x + v.x, y + v.y, z + v.z);
    }

    Point &operator+=(const Vector &v) {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    Vector operator-(const Point &p) const {
        Assert(!p.HasNaNs());
        return Vector(x - p.x, y - p.y, z - p.z);
    }

    Point operator-(const Vector &v) const {
        Assert(!v.HasNaNs());
        return Point(x - v.x, y - v.y, z - v.z);
    }

    Point &operator-=(const Vector &v) {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    Point &operator+=(const Point &p) {
        Assert(!p.HasNaNs());
        x += p.x; y += p.y; z += p.z;
        return *this;
    }
    Point operator+(const Point &p) const {
        Assert(!p.HasNaNs());
        return Point(x + p.x, y + p.y, z + p.z);
    }
    Point operator* (float f) const {
        return Point(f*x, f*y, f*z);
    }
    Point &operator*=(float f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    Point operator/ (float f) const {
        float inv = 1.f/f;
        return Point(inv*x, inv*y, inv*z);
    }
    Point &operator/=(float f) {
        float inv = 1.f/f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    float operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    float &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    bool HasNaNs() const {
        return isnan(x) || isnan(y) || isnan(z);
    }

    bool operator==(const Point &p) const {
        return x == p.x && y == p.y && z == p.z;
    }
    bool operator!=(const Point &p) const {
        return x != p.x || y != p.y || z != p.z;
    }

    // Point Public Data
    float x, y, z;
};


class Point2D;


class Vector2D {
public:
    Vector2D() : Vector2D(0, 0) {}

    explicit Vector2D(const Point2D &p);

    Vector2D(float x, float y) : x(x), y(y) { Assert(!HasNaNs()); }

#ifndef NDEBUG
    Vector2D(const Vector2D &v) : Vector2D(v.x, v.y) { Assert(!v.HasNaNs()); }

    Vector2D &operator=(const Vector2D &v) {
        Assert(!v.HasNaNs());
        x = v.x; y = v.y;
        return *this;
    }
#endif  // NDEBUG

    bool HasNaNs() const { return isnan(x) || isnan(y); }

    float LengthSquared() const { return x * x + y * y; }

    float Length() const { return sqrtf(LengthSquared()); }

    float Dot(const Vector2D &v) const {
        Assert(!v.HasNaNs());
        return x * v.x + y * v.y;
    }

    float Cross(const Vector2D &v) const {
        Assert(!v.HasNaNs());
        return x * v.y - y * v.x;
    }

    Vector2D Rotate(float angle) const {
        Vector2D x(cosf(-angle), sinf(-angle));
        return Vector2D(x.Dot(*this), x.Cross(*this));
    }

    Vector2D operator+(const Vector2D &v) const {
        Assert(!v.HasNaNs());
        return Vector2D(x + v.x, y + v.y);
    }

    Vector2D& operator+=(const Vector2D &v) {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y;
        return *this;
    }

    Vector2D operator-(const Vector2D &v) const {
        Assert(!v.HasNaNs());
        return Vector2D(x - v.x, y - v.y);
    }

    Vector2D& operator-=(const Vector2D &v) {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y;
        return *this;
    }

    Vector2D operator*(float f) const { return Vector2D(f * x, f * y); }

    Vector2D &operator*=(float f) {
        Assert(!isnan(f));
        x *= f; y *= f;
        return *this;
    }

    Vector2D operator/(float f) const {
        Assert(f != 0);
        return operator*(1.f / f);
    }

    Vector2D &operator/=(float f) {
        Assert(f != 0);
        return operator*=(1.f / f);
    }

    Vector2D operator-() const { return Vector2D(-x, -y); }

    Vector2D operator~() const { return Vector2D(-y, x); }

    float operator[](int i) const { Assert(0 <= i && i <= 1); return (&x)[i]; }

    float &operator[](int i) { Assert(0 <= i && i <= 1); return (&x)[i]; }

    bool operator==(const Vector2D &v) const { return x == v.x && y == v.y; }

    bool operator!=(const Vector2D &v) const { return x != v.x || y != v.y; }

    float x, y;
};


class Point2D {
public:
    Point2D() : Point2D(0, 0) {}

    Point2D(float x, float y) : x(x), y(y) { Assert(!HasNaNs()); }

    explicit Point2D(const Vector2D &v) : Point2D(v.x, v.y) {}

#ifndef NDEBUG
    Point2D(const Point2D &p) : Point2D(p.x, p.y) { Assert(!p.HasNaNs()); }

    Point2D &operator=(const Point2D &p) {
        Assert(!p.HasNaNs());
        x = p.x; y = p.y;
        return *this;
    }
#endif  // NDEBUG

    bool HasNaNs() const { return isnan(x) || isnan(y); }

    Point2D operator+(const Vector2D &v) const {
        Assert(!v.HasNaNs());
        return Point2D(x + v.x, y + v.y);
    }

    Point2D &operator+=(const Vector2D &v) {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y;
        return *this;
    }

    Vector2D operator-(const Point2D &p) const {
        Assert(!p.HasNaNs());
        return Vector2D(x - p.x, y - p.y);
    }

    Point2D operator-(const Vector2D &v) const {
        Assert(!v.HasNaNs());
        return Point2D(x - v.x, y - v.y);
    }

    Point2D &operator-=(const Vector2D &v) {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y;
        return *this;
    }

    Point2D &operator+=(const Point2D &p) {
        Assert(!p.HasNaNs());
        x += p.x; y += p.y;
        return *this;
    }

    Point2D operator+(const Point2D &p) const {
        Assert(!p.HasNaNs());
        return Point2D(x + p.x, y + p.y);
    }

    Point2D operator*(float f) const { return Point2D(f * x, f * y); }

    Point2D &operator*=(float f) { x *= f; y *= f; return *this; }

    Point2D operator/(float f) const { return operator*(1.f / f); }

    Point2D &operator/=(float f) { return operator*=(1.f / f); }

    float operator[](int i) const { Assert(0 <= i && i <= 1); return (&x)[i]; }

    float &operator[](int i) { Assert(0 <= i && i <= 1); return (&x)[i]; }

    bool operator==(const Point2D &p) const { return x == p.x && y == p.y; }

    bool operator!=(const Point2D &p) const { return x != p.x || y != p.y; }

    float x, y;
};


inline Vector2D::Vector2D(const Point2D &p) : Vector2D(p.x, p.y) {}


class Line2D {
public:
    Line2D() : Line2D(Point2D(0, 0), Vector2D(0, 0)) {}

    Line2D(const Point2D &p0, const Vector2D &v) : p0(p0), v(v) {}

    Point2D p0;
    Vector2D v;
};


class Circle {
public:
    Circle() : Circle(Point2D(0, 0), 0) {}

    Circle(const Point2D &center, float radius) :
            center(center), radius(radius) {}

    Point2D center;
    float radius;
};


class Arc : public Circle {
public:
    Arc() : Arc(Point2D(0, 0), 0, 0, 0) {}

    Arc(const Point2D &center, float radius, float theta0, float theta) :
            Circle(center, radius), theta0(theta0), theta(theta) {}

    bool Split(float t, Arc *a1, Arc *a2) const {
        while (t < theta0) {
            t += M_PI * 2.f;
        }
        if (t == theta0 || t >= theta0 + theta) {
            return false;
        }
        *a1 = Arc(center, radius, theta0, t - theta0);
        *a2 = Arc(center, radius, t, theta0 + theta - t);
        return true;
    }

    float theta0;
    float theta;
};


class Normal {
public:
    // Normal Public Methods
    Normal() { x = y = z = 0.f; }
    Normal(float xx, float yy, float zz)
        : x(xx), y(yy), z(zz) {
        Assert(!HasNaNs());
    }
    Normal operator-() const {
        return Normal(-x, -y, -z);
    }
    Normal operator+ (const Normal &n) const {
        Assert(!n.HasNaNs());
        return Normal(x + n.x, y + n.y, z + n.z);
    }

    Normal& operator+=(const Normal &n) {
        Assert(!n.HasNaNs());
        x += n.x; y += n.y; z += n.z;
        return *this;
    }
    Normal operator- (const Normal &n) const {
        Assert(!n.HasNaNs());
        return Normal(x - n.x, y - n.y, z - n.z);
    }

    Normal& operator-=(const Normal &n) {
        Assert(!n.HasNaNs());
        x -= n.x; y -= n.y; z -= n.z;
        return *this;
    }
    bool HasNaNs() const {
        return isnan(x) || isnan(y) || isnan(z);
    }
    Normal operator*(float f) const {
        return Normal(f*x, f*y, f*z);
    }

    Normal &operator*=(float f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    Normal operator/(float f) const {
        Assert(f != 0);
        float inv = 1.f/f;
        return Normal(x * inv, y * inv, z * inv);
    }

    Normal &operator/=(float f) {
        Assert(f != 0);
        float inv = 1.f/f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    float LengthSquared() const { return x*x + y*y + z*z; }
    float Length() const        { return sqrtf(LengthSquared()); }

#ifndef NDEBUG
    Normal(const Normal &n) {
        Assert(!n.HasNaNs());
        x = n.x; y = n.y; z = n.z;
    }

    Normal &operator=(const Normal &n) {
        Assert(!n.HasNaNs());
        x = n.x; y = n.y; z = n.z;
        return *this;
    }
#endif // !NDEBUG
    explicit Normal(const Vector &v)
      : x(v.x), y(v.y), z(v.z) {
        Assert(!v.HasNaNs());
    }
    float operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    float &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }

    bool operator==(const Normal &n) const {
        return x == n.x && y == n.y && z == n.z;
    }
    bool operator!=(const Normal &n) const {
        return x != n.x || y != n.y || z != n.z;
    }

    // Normal Public Data
    float x, y, z;
};


class Ray {
public:
    // Ray Public Methods
    Ray() : mint(0.f), maxt(INFINITY), time(0.f), depth(0) { }
    Ray(const Point &origin, const Vector &direction,
        float start, float end = INFINITY, float t = 0.f, int d = 0)
        : o(origin), d(direction), mint(start), maxt(end), time(t), depth(d) { }
    Ray(const Point &origin, const Vector &direction, const Ray &parent,
        float start, float end = INFINITY)
        : o(origin), d(direction), mint(start), maxt(end),
          time(parent.time), depth(parent.depth+1) { }
    Point operator()(float t) const { return o + d * t; }
    bool HasNaNs() const {
        return (o.HasNaNs() || d.HasNaNs() ||
                isnan(mint) || isnan(maxt));
    }

    // Ray Public Data
    Point o;
    Vector d;
    mutable float mint, maxt;
    float time;
    int depth;
};


class RayDifferential : public Ray {
public:
    // RayDifferential Public Methods
    RayDifferential() { hasDifferentials = false; }
    RayDifferential(const Point &org, const Vector &dir, float start,
        float end = INFINITY, float t = 0.f, int d = 0)
            : Ray(org, dir, start, end, t, d) {
        hasDifferentials = false;
    }
    RayDifferential(const Point &org, const Vector &dir, const Ray &parent,
        float start, float end = INFINITY)
            : Ray(org, dir, start, end, parent.time, parent.depth+1) {
        hasDifferentials = false;
    }
    explicit RayDifferential(const Ray &ray) : Ray(ray) {
        hasDifferentials = false;
    }
    bool HasNaNs() const {
        return Ray::HasNaNs() ||
           (hasDifferentials && (rxOrigin.HasNaNs() || ryOrigin.HasNaNs() ||
                                 rxDirection.HasNaNs() || ryDirection.HasNaNs()));
    }
    void ScaleDifferentials(float s) {
        rxOrigin = o + (rxOrigin - o) * s;
        ryOrigin = o + (ryOrigin - o) * s;
        rxDirection = d + (rxDirection - d) * s;
        ryDirection = d + (ryDirection - d) * s;
    }

    // RayDifferential Public Data
    bool hasDifferentials;
    Point rxOrigin, ryOrigin;
    Vector rxDirection, ryDirection;
};


class BBox {
public:
    // BBox Public Methods
    BBox() {
        pMin = Point( INFINITY,  INFINITY,  INFINITY);
        pMax = Point(-INFINITY, -INFINITY, -INFINITY);
    }
    BBox(const Point &p) : pMin(p), pMax(p) { }
    BBox(const Point &p1, const Point &p2) {
        pMin = Point(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z));
        pMax = Point(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z));
    }
    friend BBox Union(const BBox &b, const Point &p);
    friend BBox Union(const BBox &b, const BBox &b2);
    bool Overlaps(const BBox &b) const {
        bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
        bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
        bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
        return (x && y && z);
    }
    bool Inside(const Point &pt) const {
        return (pt.x >= pMin.x && pt.x <= pMax.x &&
                pt.y >= pMin.y && pt.y <= pMax.y &&
                pt.z >= pMin.z && pt.z <= pMax.z);
    }
    void Expand(float delta) {
        pMin -= Vector(delta, delta, delta);
        pMax += Vector(delta, delta, delta);
    }
    float SurfaceArea() const {
        Vector d = pMax - pMin;
        return 2.f * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    float Volume() const {
        Vector d = pMax - pMin;
        return d.x * d.y * d.z;
    }
    int MaximumExtent() const {
        Vector diag = pMax - pMin;
        if (diag.x > diag.y && diag.x > diag.z)
            return 0;
        else if (diag.y > diag.z)
            return 1;
        else
            return 2;
    }
    const Point &operator[](int i) const;
    Point &operator[](int i);
    Point Lerp(float tx, float ty, float tz) const {
        return Point(::Lerp(tx, pMin.x, pMax.x), ::Lerp(ty, pMin.y, pMax.y),
                     ::Lerp(tz, pMin.z, pMax.z));
    }
    Vector Offset(const Point &p) const {
        return Vector((p.x - pMin.x) / (pMax.x - pMin.x),
                      (p.y - pMin.y) / (pMax.y - pMin.y),
                      (p.z - pMin.z) / (pMax.z - pMin.z));
    }
    void BoundingSphere(Point *c, float *rad) const;
    bool IntersectP(const Ray &ray, float *hitt0 = NULL, float *hitt1 = NULL) const;

    bool operator==(const BBox &b) const {
        return b.pMin == pMin && b.pMax == pMax;
    }
    bool operator!=(const BBox &b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }

    // BBox Public Data
    Point pMin, pMax;
};



// Geometry Inline Functions
inline Vector::Vector(const Point &p)
    : x(p.x), y(p.y), z(p.z) {
    Assert(!HasNaNs());
}


inline Vector operator*(float f, const Vector &v) { return v*f; }

inline Vector2D operator*(float f, const Vector2D &v) { return v * f; }

inline float Dot(const Vector &v1, const Vector &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline float Dot(const Vector2D &v1, const Vector2D &v2) {
    return v1.Dot(v2);
}


inline float AbsDot(const Vector &v1, const Vector &v2) {
    return fabsf(Dot(v1, v2));
}


inline Vector Cross(const Vector &v1, const Vector &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


inline float Cross(const Vector2D &v1, const Vector2D &v2) {
    return v1.Cross(v2);
}


inline Vector Cross(const Vector &v1, const Normal &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


inline Vector Cross(const Normal &v1, const Vector &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


inline Vector Normalize(const Vector &v) { return v / v.Length(); }

inline Vector2D Normalize(const Vector2D &v) { return v / v.Length(); }

inline void CoordinateSystem(const Vector &v1, Vector *v2, Vector *v3) {
    if (fabsf(v1.x) > fabsf(v1.y)) {
        float invLen = 1.f / sqrtf(v1.x*v1.x + v1.z*v1.z);
        *v2 = Vector(-v1.z * invLen, 0.f, v1.x * invLen);
    }
    else {
        float invLen = 1.f / sqrtf(v1.y*v1.y + v1.z*v1.z);
        *v2 = Vector(0.f, v1.z * invLen, -v1.y * invLen);
    }
    *v3 = Cross(v1, *v2);
}


inline float Distance(const Point &p1, const Point &p2) {
    return (p1 - p2).Length();
}

inline float Distance(const Point2D &p1, const Point2D &p2) {
    return (p1 - p2).Length();
}


inline float DistanceSquared(const Point &p1, const Point &p2) {
    return (p1 - p2).LengthSquared();
}


inline float DistanceSquared(const Point2D &p1, const Point2D &p2) {
    return (p1 - p2).LengthSquared();
}


inline Point operator*(float f, const Point &p) {
    Assert(!p.HasNaNs());
    return p*f;
}

inline Point2D operator*(float f, const Point2D &p) {
    Assert(!p.HasNaNs());
    return p * f;
}


inline Normal operator*(float f, const Normal &n) {
    return Normal(f*n.x, f*n.y, f*n.z);
}


inline Normal Normalize(const Normal &n) {
    return n / n.Length();
}


inline Vector::Vector(const Normal &n)
  : x(n.x), y(n.y), z(n.z) {
    Assert(!n.HasNaNs());
}


inline float Dot(const Normal &n1, const Vector &v2) {
    Assert(!n1.HasNaNs() && !v2.HasNaNs());
    return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}


inline float Dot(const Vector &v1, const Normal &n2) {
    Assert(!v1.HasNaNs() && !n2.HasNaNs());
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}


inline float Dot(const Normal &n1, const Normal &n2) {
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}


inline float AbsDot(const Normal &n1, const Vector &v2) {
    Assert(!n1.HasNaNs() && !v2.HasNaNs());
    return fabsf(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}


inline float AbsDot(const Vector &v1, const Normal &n2) {
    Assert(!v1.HasNaNs() && !n2.HasNaNs());
    return fabsf(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}


inline float AbsDot(const Normal &n1, const Normal &n2) {
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
    return fabsf(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}


inline Normal Faceforward(const Normal &n, const Vector &v) {
    return (Dot(n, v) < 0.f) ? -n : n;
}


inline Normal Faceforward(const Normal &n, const Normal &n2) {
    return (Dot(n, n2) < 0.f) ? -n : n;
}



inline Vector Faceforward(const Vector &v, const Vector &v2) {
    return (Dot(v, v2) < 0.f) ? -v : v;
}



inline Vector Faceforward(const Vector &v, const Normal &n2) {
    return (Dot(v, n2) < 0.f) ? -v : v;
}


inline const Point &BBox::operator[](int i) const {
    Assert(i == 0 || i == 1);
    return (&pMin)[i];
}



inline Point &BBox::operator[](int i) {
    Assert(i == 0 || i == 1);
    return (&pMin)[i];
}


inline Vector SphericalDirection(float sintheta,
                                 float costheta, float phi) {
    return Vector(sintheta * cosf(phi),
                  sintheta * sinf(phi),
                  costheta);
}


inline Vector SphericalDirection(float sintheta, float costheta,
                                 float phi, const Vector &x,
                                 const Vector &y, const Vector &z) {
    return sintheta * cosf(phi) * x +
           sintheta * sinf(phi) * y + costheta * z;
}


inline float SphericalTheta(const Vector &v) {
    return acosf(Clamp(v.z, -1.f, 1.f));
}


inline float SphericalPhi(const Vector &v) {
    float p = atan2f(v.y, v.x);
    return (p < 0.f) ? p + 2.f*M_PI : p;
}



#endif // PBRT_CORE_GEOMETRY_H
