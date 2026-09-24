// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef STL_READER_H_
#define STL_READER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <iomanip>   // for setprecision

// ============================================================
// STL geometry: reading a triangle mesh, placing it, and answering
// point-in-solid and signed-distance queries with a BVH.
//
// A case uses it by including this header in its Parm.H, reading a file with
// ReadSTLFile, placing the mesh with CenterAt / NormalizeSize / RotateAroundPoint
// followed by ApplyTransforms, calling BuildBVH once, and then querying with
// IsPointInsideObject and SignedDistanceBVH. See exec/_Tests/COMPAS-STL-5Eq.
//
// Host only. Nothing here carries an AMREX_GPU_ annotation, the queries allocate
// a traversal stack, and the mesh lives in host memory, so a case that calls
// these from initial_condition must be built for the CPU.
//
// Geometry is stored in float, and the queries assume an O(1) mesh, which is what
// NormalizeSize gives. Distances are returned positive inside the solid.
// ============================================================

namespace StlGeom {

// ============================================================
// Basic Data Structures
// ============================================================

struct Vector3 {
    float x, y, z;
};
static_assert(sizeof(Vector3) == 12, "Vector3 must be tightly packed (3 floats).");

struct Triangle {
    Vector3 vertices[3];
};

// M_PI is POSIX rather than standard C++, so keep a local constant
constexpr float PI_F = 3.14159265358979323846f;

// ============================================================
// STLObject
// ============================================================

struct STLObject {

    std::string filename;
    std::vector<Triangle> triangles;

    // ============================================================
    // Read an STL file: binary first, ASCII if the header says so, then a
    // best-effort binary pass for files whose size does not match the count.
    // Returns false and leaves the mesh empty on any failure.
    // ============================================================
    bool ReadSTLFile()
    {
        std::ifstream stl(filename, std::ios::binary);
        if (!stl) {
            std::cerr << "Cannot open " << filename << "\n";
            return false;
        }

        // Get file size for validation
        stl.seekg(0, std::ios::end);
        std::streamoff filesize = stl.tellg();
        stl.seekg(0, std::ios::beg);

        // Read header
        char header[80];
        stl.read(header, 80);
        if (!stl) { std::cerr << "Failed to read STL header\n"; return false; }

        // Try reading binary count (uint32)
        uint32_t N = 0;
        stl.read(reinterpret_cast<char*>(&N), sizeof(N));
        if (!stl) { std::cerr << "Failed to read triangle count\n"; return false; }

        // Validate expected binary size: 80 + 4 + N*50
        const std::streamoff expectedBinSize = 80 + 4 + static_cast<std::streamoff>(N) * 50;
        if (filesize == expectedBinSize) {
            if (N == 0) { std::cerr << filename << " is a binary STL with no triangles\n"; return false; }
            // Likely a binary STL, parse it safely
            triangles.clear();
            triangles.resize(N);

            for (uint32_t t = 0; t < N; ++t)
            {
                float normal[3];
                stl.read(reinterpret_cast<char*>(normal), sizeof(normal));
                if (!stl) { std::cerr << "Unexpected EOF reading normal at triangle " << t << "\n"; return false; }

                // Read 3 vertices
                for (int v = 0; v < 3; ++v) {
                    stl.read(reinterpret_cast<char*>(&triangles[t].vertices[v]), sizeof(Vector3));
                    if (!stl) { std::cerr << "Unexpected EOF reading vertex at triangle " << t << " vertex " << v << "\n"; return false; }
                }

                // attribute byte count (2 bytes), ignored, but the stream must advance
                uint16_t attr = 0;
                stl.read(reinterpret_cast<char*>(&attr), sizeof(attr));
                if (!stl) { std::cerr << "Unexpected EOF reading attribute bytes at triangle " << t << "\n"; return false; }
            }
            return true;
        }

        // If we get here, file size didn't match expected binary size.
        // Heuristic: if header starts with "solid" -> try ASCII parser fallback.
        std::string hdrStr(header, header + 80);
        auto trim = [](const std::string &s)->std::string {
            size_t a = s.find_first_not_of(" \t\r\n");
            size_t b = s.find_last_not_of(" \t\r\n");
            return (a==std::string::npos) ? std::string() : s.substr(a, b-a+1);
        };
        if (trim(hdrStr).rfind("solid", 0) == 0) {
            // Attempt simple ASCII STL parse
            stl.close();
            std::ifstream txt(filename);
            if (!txt) { std::cerr << "Cannot re-open file for ASCII parse\n"; return false; }

            triangles.clear();
            Triangle cur;
            int vcount = 0;
            std::string token;
            while (txt >> token) {
                if (token == "vertex") {
                    float x,y,z;
                    if (!(txt >> x >> y >> z)) {
                        std::cerr << "Malformed vertex in ASCII STL\n";
                        return false;
                    }
                    cur.vertices[vcount++] = {x,y,z};
                    if (vcount == 3) {
                        triangles.push_back(cur);
                        vcount = 0;
                    }
                }
                // ignore other tokens
            }
            if (triangles.empty()) {
                std::cerr << "ASCII STL parse produced no triangles\n";
                return false;
            }
            return true;
        }

        // If not ASCII (or ASCII heuristic failed) but size didn't match, try best-effort binary parse:
        // Rewind and parse binary using file size guards to avoid overruns.
        stl.clear();
        stl.seekg(84, std::ios::beg); // after header+count

        // N comes from the file and did not agree with its size above, so bound it
        // by the bytes that are actually there before allocating anything.
        const std::streamoff avail = filesize - 84;
        const uint64_t maxTris = (avail > 0) ? static_cast<uint64_t>(avail) / 50u : 0u;
        if (static_cast<uint64_t>(N) > maxTris) {
            std::cerr << "Binary STL claims " << N << " triangles but holds room for "
                      << maxTris << ": " << filename << "\n";
            return false;
        }
        triangles.clear();
        triangles.resize(N);

        for (uint32_t t = 0; t < N; ++t)
        {
            // convert tellg() result to numeric offset
            std::streampos posp = stl.tellg();
            if (posp == std::streampos(-1)) {
                std::cerr << "tellg() failed while parsing binary STL at triangle " << t << "\n";
                triangles.clear();
                return false;
            }
            std::streamoff pos = static_cast<std::streamoff>(posp);
            if (pos + static_cast<std::streamoff>(50) > filesize) {
                std::cerr << "Binary STL shorter than expected. Triangles read so far: " << t << "\n";
                triangles.clear();
                return false;
            }

            float normal[3];
            stl.read(reinterpret_cast<char*>(normal), sizeof(normal));
            for (int v = 0; v < 3; ++v) {
                stl.read(reinterpret_cast<char*>(&triangles[t].vertices[v]), sizeof(Vector3));
            }
            uint16_t attr;
            stl.read(reinterpret_cast<char*>(&attr), sizeof(attr));
            if (!stl) {
                std::cerr << "Read failed at triangle " << t << " of " << filename << "\n";
                triangles.clear();
                return false;
            }
        }
        if (triangles.empty()) {
            std::cerr << "No triangles read from " << filename << "\n";
            return false;
        }
        return true;
    }

    // ============================================================
    // 4x4 transform matrix system (column-major)
    // ============================================================

    float M[16];   // main transform matrix (column-major)

    static void MatIdentity(float A[16]) {
        for (int i = 0; i < 16; ++i) A[i] = 0.0f;
        A[0] = A[5] = A[10] = A[15] = 1.0f;
    }

    static void MatMul(float R[16], const float A[16], const float B[16]) {
        // Column-major multiplication: R = A * B
        float T[16];
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                // element at row r, col c stored at index c*4 + r
                T[c*4 + r] =
                    A[0*4 + r] * B[c*4 + 0] +
                    A[1*4 + r] * B[c*4 + 1] +
                    A[2*4 + r] * B[c*4 + 2] +
                    A[3*4 + r] * B[c*4 + 3];
            }
        }
        std::memcpy(R, T, 16 * sizeof(float));
    }

    static Vector3 TransformPoint(const float mat[16], const Vector3& v) {
        // column-major: x' = m00*x + m01*y + m02*z + m03
        return {
            mat[0]*v.x + mat[4]*v.y + mat[8]*v.z  + mat[12],
            mat[1]*v.x + mat[5]*v.y + mat[9]*v.z  + mat[13],
            mat[2]*v.x + mat[6]*v.y + mat[10]*v.z + mat[14]
        };
    }

    // -----------------------
    // Reset transform
    // -----------------------
    void ResetTransform() { MatIdentity(M); }

    // -----------------------
    // Translate, Scale, Rotate
    // (Each operation pre-multiplies the current matrix:
    //  M := Op * M)
    // -----------------------
    void Translate(Vector3 t) {
        float T[16]; MatIdentity(T);
        T[12] = t.x; T[13] = t.y; T[14] = t.z;
        MatMul(M, T, M);
    }

    void Scale(Vector3 s) {
        float S[16]; MatIdentity(S);
        S[0] = s.x; S[5] = s.y; S[10] = s.z;
        MatMul(M, S, M);
    }

    void UniformScale(float s) { Scale({s, s, s}); }

    void RotateX(float deg) {
        float a = deg * PI_F / 180.0f, c = cosf(a), s = sinf(a);
        float R[16]; MatIdentity(R);
        R[5] = c; R[9] = -s;
        R[6] = s; R[10] = c;
        MatMul(M, R, M);
    }
    void RotateY(float deg) {
        float a = deg * PI_F / 180.0f, c = cosf(a), s = sinf(a);
        float R[16]; MatIdentity(R);
        R[0] = c; R[8] = s;
        R[2] = -s; R[10] = c;
        MatMul(M, R, M);
    }
    void RotateZ(float deg) {
        float a = deg * PI_F / 180.0f, c = cosf(a), s = sinf(a);
        float R[16]; MatIdentity(R);
        R[0] = c; R[4] = -s;
        R[1] = s; R[5] = c;
        MatMul(M, R, M);
    }
    void RotateXYZ(Vector3 d) {
        RotateX(d.x); RotateY(d.y); RotateZ(d.z);
    }

    // ============================================================
    // Rotate Around Arbitrary Point
    // ============================================================
    void RotateAroundPoint(Vector3 P, Vector3 deg)
    {
        Translate({-P.x, -P.y, -P.z});
        RotateXYZ(deg);
        Translate(P);
    }

    // ============================================================
    // Internal bounding box
    // ============================================================
    void ComputeBounds(Vector3& mn, Vector3& mx) const
    {
        if (triangles.empty()) {
            mn = {0.f, 0.f, 0.f};
            mx = {0.f, 0.f, 0.f};
            return;
        }

        mn = {  std::numeric_limits<float>::infinity(),
                std::numeric_limits<float>::infinity(),
                std::numeric_limits<float>::infinity() };
        mx = { -std::numeric_limits<float>::infinity(),
               -std::numeric_limits<float>::infinity(),
               -std::numeric_limits<float>::infinity() };

        for (const auto& T : triangles) {
            for (int i = 0; i < 3; ++i) {
                const Vector3& v = T.vertices[i];
                mn.x = std::min(mn.x, v.x);
                mn.y = std::min(mn.y, v.y);
                mn.z = std::min(mn.z, v.z);

                mx.x = std::max(mx.x, v.x);
                mx.y = std::max(mx.y, v.y);
                mx.z = std::max(mx.z, v.z);
            }
        }
    }

    // ============================================================
    // Fit Into Box
    // ============================================================
    void FitIntoBox(Vector3 mn_target, Vector3 mx_target)
    {
        Vector3 mn, mx;
        ComputeBounds(mn, mx);

        Vector3 size = { mx.x - mn.x, mx.y - mn.y, mx.z - mn.z };
        // protect against zero-size components
        const float EPS = 1e-9f;
        if (size.x < EPS) size.x = EPS;
        if (size.y < EPS) size.y = EPS;
        if (size.z < EPS) size.z = EPS;

        Vector3 tgt  = { mx_target.x - mn_target.x,
                         mx_target.y - mn_target.y,
                         mx_target.z - mn_target.z };

        float sx = tgt.x / size.x;
        float sy = tgt.y / size.y;
        float sz = tgt.z / size.z;
        float s  = std::min(std::min(sx, sy), sz); // uniform

        // Center object into target box
        Vector3 oldC = {0.5f*(mn.x + mx.x), 0.5f*(mn.y + mx.y), 0.5f*(mn.z + mx.z)};
        Vector3 newC = {0.5f*(mn_target.x + mx_target.x),
                        0.5f*(mn_target.y + mx_target.y),
                        0.5f*(mn_target.z + mx_target.z)};

        Translate({-oldC.x,-oldC.y,-oldC.z});
        UniformScale(s);
        Translate(newC);
    }

    // ============================================================
    // Normalize Size (max dimension)
    // ============================================================
    void NormalizeSize(float maxSize=1.0f)
    {
        Vector3 mn, mx;
        ComputeBounds(mn, mx);
        float sx = mx.x - mn.x;
        float sy = mx.y - mn.y;
        float sz = mx.z - mn.z;
        const float EPS = 1e-9f;
        float L = std::max({sx, sy, sz, EPS});
        float s = maxSize / L;
        UniformScale(s);
    }

    // ============================================================
    // Center at a given point
    // ============================================================
    void CenterAt(Vector3 C)
    {
        Vector3 mn, mx;
        ComputeBounds(mn, mx);
        Vector3 oldC = {0.5f*(mn.x + mx.x), 0.5f*(mn.y + mx.y), 0.5f*(mn.z + mx.z)};
        Translate({C.x - oldC.x, C.y - oldC.y, C.z - oldC.z});
    }

    // ============================================================
    // Bake matrix into geometry
    // ============================================================
    void ApplyTransforms()
    {
        for (auto& t : triangles) {
            for (int i = 0; i < 3; ++i)
                t.vertices[i] = TransformPoint(M, t.vertices[i]);
        }
        ResetTransform();
    }

    // ============================================================
    // BVH acceleration structure
    // ============================================================

    struct AABB {
        Vector3 mn, mx;
    };

    struct BVHNode {
        AABB box;
        int left  = -1;
        int right = -1;
        int start = 0;   // index of first triangle in leaf
        int count = 0;   // number of triangles (0 if internal node)
    };

    std::vector<BVHNode> bvh;
    std::vector<int> primIndices;   // triangle index list (permutation)

    // ------------------------------------------------------------
    // AABB utilities
    // ------------------------------------------------------------
    static AABB TriBounds(const Triangle& T)
    {
        AABB a;
        a.mn = { std::numeric_limits<float>::infinity(),
                 std::numeric_limits<float>::infinity(),
                 std::numeric_limits<float>::infinity() };
        a.mx = { -std::numeric_limits<float>::infinity(),
                 -std::numeric_limits<float>::infinity(),
                 -std::numeric_limits<float>::infinity() };

        for (int i = 0; i < 3; ++i) {
            const Vector3& v = T.vertices[i];
            a.mn.x = std::min(a.mn.x, v.x);
            a.mn.y = std::min(a.mn.y, v.y);
            a.mn.z = std::min(a.mn.z, v.z);

            a.mx.x = std::max(a.mx.x, v.x);
            a.mx.y = std::max(a.mx.y, v.y);
            a.mx.z = std::max(a.mx.z, v.z);
        }
        return a;
    }

    static AABB Merge(const AABB& a, const AABB& b)
    {
        AABB r;
        r.mn.x = std::min(a.mn.x, b.mn.x);
        r.mn.y = std::min(a.mn.y, b.mn.y);
        r.mn.z = std::min(a.mn.z, b.mn.z);

        r.mx.x = std::max(a.mx.x, b.mx.x);
        r.mx.y = std::max(a.mx.y, b.mx.y);
        r.mx.z = std::max(a.mx.z, b.mx.z);

        return r;
    }

    // ------------------------------------------------------------
    // Robust ray-AABB test (slab method), safe for dir components near zero
    // ------------------------------------------------------------
    static inline bool RayAABB(const Vector3& o, const Vector3& d,
                               const AABB& B)
    {
        const float EPS = 1e-9f;
        float tmin = -std::numeric_limits<float>::infinity();
        float tmax =  std::numeric_limits<float>::infinity();

        // X
        if (std::fabs(d.x) < EPS) {
            if (o.x < B.mn.x || o.x > B.mx.x) return false;
        } else {
            float tx1 = (B.mn.x - o.x) / d.x;
            float tx2 = (B.mx.x - o.x) / d.x;
            if (tx1 > tx2) std::swap(tx1, tx2);
            tmin = std::max(tmin, tx1);
            tmax = std::min(tmax, tx2);
            if (tmin > tmax) return false;
        }

        // Y
        if (std::fabs(d.y) < EPS) {
            if (o.y < B.mn.y || o.y > B.mx.y) return false;
        } else {
            float ty1 = (B.mn.y - o.y) / d.y;
            float ty2 = (B.mx.y - o.y) / d.y;
            if (ty1 > ty2) std::swap(ty1, ty2);
            tmin = std::max(tmin, ty1);
            tmax = std::min(tmax, ty2);
            if (tmin > tmax) return false;
        }

        // Z
        if (std::fabs(d.z) < EPS) {
            if (o.z < B.mn.z || o.z > B.mx.z) return false;
        } else {
            float tz1 = (B.mn.z - o.z) / d.z;
            float tz2 = (B.mx.z - o.z) / d.z;
            if (tz1 > tz2) std::swap(tz1, tz2);
            tmin = std::max(tmin, tz1);
            tmax = std::min(tmax, tz2);
            if (tmin > tmax) return false;
        }

        return true;
    }

    // ------------------------------------------------------------
    // Moller-Trumbore ray-triangle intersection (writes t to outT if outT != nullptr)
    // ------------------------------------------------------------
    static inline bool RayTriIntersect(const Vector3& orig,
                                       const Vector3& dir,
                                       const Triangle& tri,
                                       float* outT = nullptr,
                                       float* outU = nullptr,
                                       float* outV = nullptr)
    {
        const Vector3& v0 = tri.vertices[0];
        const Vector3& v1 = tri.vertices[1];
        const Vector3& v2 = tri.vertices[2];

        Vector3 e1{ v1.x - v0.x, v1.y - v0.y, v1.z - v0.z };
        Vector3 e2{ v2.x - v0.x, v2.y - v0.y, v2.z - v0.z };

        // p = cross(dir, e2)
        Vector3 p{
            dir.y * e2.z - dir.z * e2.y,
            dir.z * e2.x - dir.x * e2.z,
            dir.x * e2.y - dir.y * e2.x
        };

        float det = e1.x * p.x + e1.y * p.y + e1.z * p.z;

        // det scales as (edge length)^2 times |dir|, so testing it against a fixed
        // absolute epsilon would reject every triangle of a small mesh. Compare it
        // instead with the magnitudes it is built from, which is dimensionless.
        const float e1e1 = e1.x*e1.x + e1.y*e1.y + e1.z*e1.z;
        const float pp   = p.x*p.x + p.y*p.y + p.z*p.z;
        if (det * det <= 1e-12f * e1e1 * pp) return false;   // ray parallel to the triangle

        float invDet = 1.0f / det;

        Vector3 s{
            orig.x - v0.x,
            orig.y - v0.y,
            orig.z - v0.z
        };

        float u = invDet * (s.x * p.x + s.y * p.y + s.z * p.z);
        if (u < 0.0f || u > 1.0f) return false;

        // q = cross(s, e1)
        Vector3 q{
            s.y * e1.z - s.z * e1.y,
            s.z * e1.x - s.x * e1.z,
            s.x * e1.y - s.y * e1.x
        };

        float v = invDet * (dir.x * q.x + dir.y * q.y + dir.z * q.z);
        if (v < 0.0f || u + v > 1.0f) return false;

        float t = invDet * (e2.x * q.x + e2.y * q.y + e2.z * q.z);
        if (outT) *outT = t;
        if (outU) *outU = u;
        if (outV) *outV = v;
        return t > 0.0f;
    }

    // ------------------------------------------------------------
    // Recursive BVH Builder
    // ------------------------------------------------------------
    int BuildNode(int start, int end)
    {
        BVHNode node;
        node.start = start;
        node.count = end - start;

        // Compute bounding box of this range
        AABB box = TriBounds(triangles[primIndices[start]]);
        for (int i = start + 1; i < end; ++i)
            box = Merge(box, TriBounds(triangles[primIndices[i]]));

        node.box = box;

        int n = end - start;
        int index = static_cast<int>(bvh.size());
        bvh.push_back(node);

        // Leaf
        const int maxLeafSize = 4;
        if (n <= maxLeafSize) {
            // leaf node already stored with start/count
            return index;
        }

        // Split on longest axis (robust)
        Vector3 extent {
            box.mx.x - box.mn.x,
            box.mx.y - box.mn.y,
            box.mx.z - box.mn.z
        };

        int axis = 0;
        if (extent.y > extent.x) axis = 1;
        if (extent.z > ((axis == 1) ? extent.y : extent.x)) axis = 2;

        float mid = 0.5f * ( (axis == 0) ? (box.mn.x + box.mx.x)
                            : (axis == 1) ? (box.mn.y + box.mx.y)
                                          : (box.mn.z + box.mx.z) );

        // partition by centroid along axis
        auto midIt = std::partition(
            primIndices.begin() + start,
            primIndices.begin() + end,
            [&](int idx)
            {
                const Triangle& T = triangles[idx];
                float c;
                if (axis == 0)
                    c = (T.vertices[0].x + T.vertices[1].x + T.vertices[2].x) / 3.0f;
                else if (axis == 1)
                    c = (T.vertices[0].y + T.vertices[1].y + T.vertices[2].y) / 3.0f;
                else
                    c = (T.vertices[0].z + T.vertices[1].z + T.vertices[2].z) / 3.0f;
                return c < mid;
            }
        );

        int midIdx = static_cast<int>(std::distance(primIndices.begin(), midIt));
        // If partition failed to separate, split evenly
        if (midIdx == start || midIdx == end) midIdx = start + (n / 2);

        int L = BuildNode(start, midIdx);
        int R = BuildNode(midIdx, end);

        bvh[index].left = L;
        bvh[index].right = R;
        bvh[index].count = 0; // internal

        return index;
    }

    // ------------------------------------------------------------
    // BuildBVH(): public
    // ------------------------------------------------------------
    void BuildBVH()
    {
        primIndices.resize(triangles.size());
        for (size_t i = 0; i < triangles.size(); ++i)
            primIndices[i] = static_cast<int>(i);

        bvh.clear();
        bvh.reserve(triangles.size() * 2 + 1);

        if (!primIndices.empty())
            BuildNode(0, static_cast<int>(primIndices.size()));
    }

    // ------------------------------------------------------------
    // BVH Ray traversal (iterative with dynamic stack)
    // Returns number of triangle hits (useful for odd-even point-in-poly)
    // ------------------------------------------------------------
    int RayIntersectBVH(const Vector3& orig, const Vector3& dir) const
    {
        if (bvh.empty()) return 0;
        int hits = 0;

        // Small offset to avoid starting exactly on geometry (helps edge cases),
        // scaled to the mesh so that it survives at any size.
        const Vector3 ext { bvh[0].box.mx.x - bvh[0].box.mn.x,
                            bvh[0].box.mx.y - bvh[0].box.mn.y,
                            bvh[0].box.mx.z - bvh[0].box.mn.z };
        const float meshScale = std::max(ext.x, std::max(ext.y, ext.z));
        const float ORIGIN_EPS = 1e-6f * ((meshScale > 0.0f) ? meshScale : 1.0f);
        Vector3 origin2 = { orig.x + dir.x * ORIGIN_EPS,
                            orig.y + dir.y * ORIGIN_EPS,
                            orig.z + dir.z * ORIGIN_EPS };

        std::vector<int> stack;
        stack.reserve(64);
        stack.push_back(0); // root index

        while (!stack.empty())
        {
            int idx = stack.back();
            stack.pop_back();
            if (idx < 0) continue;

            const BVHNode& N = bvh[idx];

            if (!RayAABB(origin2, dir, N.box))
                continue;

            if (N.count > 0) {
                // Leaf: test all triangles in leaf
                for (int i = 0; i < N.count; ++i) {
                    int triIdx = primIndices[N.start + i];
                    float t;
                    if (RayTriIntersect(origin2, dir, triangles[triIdx], &t)) {
                        if (t > 0.0f) ++hits;
                    }
                }
            } else {
                // Internal: push children (order doesn't matter for count)
                if (N.left >= 0)  stack.push_back(N.left);
                if (N.right >= 0) stack.push_back(N.right);
            }
        }

        return hits;
    }

    // ------------------------------------------------------------
    // Point-inside: odd-even rule (shoot ray along non-axis vector)
    // ------------------------------------------------------------
    bool IsPointInsideObject(Vector3 P) const
    {
        // Deterministic non-axis ray direction (less likely to hit edges exactly).
        Vector3 dir{ 1.234567f, 0.54321f, 0.123456f };

        // Normalize dir
        float invLen = 1.0f / std::sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
        dir.x *= invLen; dir.y *= invLen; dir.z *= invLen;

        return (RayIntersectBVH(P, dir) & 1) != 0;
    }

    // ------------------------------------------------------------
    // Diagnostics
    // ------------------------------------------------------------

    // Print basic mesh summary: count, bounds, center, size
    void PrintMeshSummary(std::ostream& os = std::cout) const
    {
        Vector3 mn, mx;
        ComputeBounds(mn, mx);
        Vector3 size = { mx.x - mn.x, mx.y - mn.y, mx.z - mn.z };
        Vector3 center = { 0.5f*(mn.x + mx.x), 0.5f*(mn.y + mx.y), 0.5f*(mn.z + mx.z) };

        os << std::fixed << std::setprecision(6);
        os << "Mesh summary:\n";
        os << "  triangles: " << triangles.size() << "\n";
        os << "  bounds min: (" << mn.x << ", " << mn.y << ", " << mn.z << ")\n";
        os << "  bounds max: (" << mx.x << ", " << mx.y << ", " << mx.z << ")\n";
        os << "  center:     (" << center.x << ", " << center.y << ", " << center.z << ")\n";
        os << "  size:       (" << size.x << ", " << size.y << ", " << size.z << ")\n";
    }

    // ------------------------- Closest-point / distance utilities -------------------------

    // Compute squared distance from point p to triangle tri, and write barycentric coords (u,v) and tParam if wanted.
    // Returns squared distance; also outputs closest point via outClosest if non-null.
    static inline float PointTriangleDistanceSq(const Vector3& p, const Triangle& tri,
                                                Vector3* outClosest = nullptr,
                                                float* outU = nullptr, float* outV = nullptr)
    {
        // from "Real-Time Collision Detection" (Christer Ericson)
        const Vector3& A = tri.vertices[0];
        const Vector3& B = tri.vertices[1];
        const Vector3& C = tri.vertices[2];

        // AB, AC, AP
        float ABx = B.x - A.x, ABy = B.y - A.y, ABz = B.z - A.z;
        float ACx = C.x - A.x, ACy = C.y - A.y, ACz = C.z - A.z;
        float APx = p.x - A.x, APy = p.y - A.y, APz = p.z - A.z;

        float d1 = ABx*APx + ABy*APy + ABz*APz;
        float d2 = ACx*APx + ACy*APy + ACz*APz;
        if (d1 <= 0.0f && d2 <= 0.0f) {
            if (outClosest) *outClosest = A;
            if (outU) *outU = 0.0f; if (outV) *outV = 0.0f;
            float dx = p.x - A.x, dy = p.y - A.y, dz = p.z - A.z;
            return dx*dx + dy*dy + dz*dz;
        }

        // Check vertex region B
        float BPx = p.x - B.x, BPy = p.y - B.y, BPz = p.z - B.z;
        float d3 = ABx*BPx + ABy*BPy + ABz*BPz;
        float d4 = ACx*BPx + ACy*BPy + ACz*BPz;
        if (d3 >= 0.0f && d4 <= d3) {
            if (outClosest) *outClosest = B;
            if (outU) *outU = 1.0f; if (outV) *outV = 0.0f;
            float dx = p.x - B.x, dy = p.y - B.y, dz = p.z - B.z;
            return dx*dx + dy*dy + dz*dz;
        }

        // Edge region AB
        float vc = d1*d4 - d3*d2;
        if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
            float v = d1 / (d1 - d3);
            if (outClosest) {
                outClosest->x = A.x + v * ABx;
                outClosest->y = A.y + v * ABy;
                outClosest->z = A.z + v * ABz;
            }
            if (outU) *outU = v; if (outV) *outV = 0.0f;
            float dx = p.x - (A.x + v * ABx);
            float dy = p.y - (A.y + v * ABy);
            float dz = p.z - (A.z + v * ABz);
            return dx*dx + dy*dy + dz*dz;
        }

        // Vertex region C
        float CPx = p.x - C.x, CPy = p.y - C.y, CPz = p.z - C.z;
        float d5 = ABx*CPx + ABy*CPy + ABz*CPz;
        float d6 = ACx*CPx + ACy*CPy + ACz*CPz;
        if (d6 >= 0.0f && d5 <= d6) {
            if (outClosest) *outClosest = C;
            if (outU) *outU = 0.0f; if (outV) *outV = 1.0f;
            float dx = p.x - C.x, dy = p.y - C.y, dz = p.z - C.z;
            return dx*dx + dy*dy + dz*dz;
        }

        // Edge region AC
        float vb = d5*d2 - d1*d6;
        if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
            float w = d2 / (d2 - d6);
            if (outClosest) {
                outClosest->x = A.x + w * ACx;
                outClosest->y = A.y + w * ACy;
                outClosest->z = A.z + w * ACz;
            }
            if (outU) *outU = 0.0f; if (outV) *outV = w;
            float dx = p.x - (A.x + w * ACx);
            float dy = p.y - (A.y + w * ACy);
            float dz = p.z - (A.z + w * ACz);
            return dx*dx + dy*dy + dz*dz;
        }

        // Edge region BC
        float va = d3*d6 - d5*d4;
        if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
            float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            if (outClosest) {
                outClosest->x = B.x + w * (C.x - B.x);
                outClosest->y = B.y + w * (C.y - B.y);
                outClosest->z = B.z + w * (C.z - B.z);
            }
            if (outU) *outU = 1.0f - w; if (outV) *outV = w;
            float cx = B.x + w * (C.x - B.x);
            float cy = B.y + w * (C.y - B.y);
            float cz = B.z + w * (C.z - B.z);
            float dx = p.x - cx, dy = p.y - cy, dz = p.z - cz;
            return dx*dx + dy*dy + dz*dz;
        }

        // Inside face region. Compute barycentric coordinates (u,v,w).
        // va + vb + vc is twice the triangle area (Ericson, Real-Time Collision
        // Detection, 5.1.5), so it vanishes for a degenerate triangle. Fall back to
        // the nearest vertex there rather than dividing by zero.
        if (!(va + vb + vc > 0.0f)) {
            float best = std::numeric_limits<float>::infinity();
            const Vector3* bestV = &A;
            for (const Vector3* V : { &A, &B, &C }) {
                float ddx = p.x - V->x, ddy = p.y - V->y, ddz = p.z - V->z;
                float dd = ddx*ddx + ddy*ddy + ddz*ddz;
                if (dd < best) { best = dd; bestV = V; }
            }
            if (outClosest) *outClosest = *bestV;
            if (outU) *outU = 0.0f;
            if (outV) *outV = 0.0f;
            return best;
        }
        float denom = 1.0f / (va + vb + vc);
        float v = vb * denom;
        float w = vc * denom;
        float u = 1.0f - v - w;
        if (outClosest) {
            outClosest->x = u*A.x + v*B.x + w*C.x;
            outClosest->y = u*A.y + v*B.y + w*C.y;
            outClosest->z = u*A.z + v*B.z + w*C.z;
        }
        if (outU) *outU = v;
        if (outV) *outV = w;
        float dx = p.x - (u*A.x + v*B.x + w*C.x);
        float dy = p.y - (u*A.y + v*B.y + w*C.y);
        float dz = p.z - (u*A.z + v*B.z + w*C.z);
        return dx*dx + dy*dy + dz*dz;
    }

    // -------------------- Nearest distance, tested against every triangle. Used as
    // the fallback when BuildBVH has not been called; the BVH version below is what
    // the queries use once it has. --------------------
    float UnsignedDistanceBrute(const Vector3& p) const
    {
        float best = std::numeric_limits<float>::infinity();
        Vector3 tmp;
        for (size_t i = 0; i < triangles.size(); ++i) {
            float dsq = PointTriangleDistanceSq(p, triangles[i], &tmp);
            if (dsq < best) best = dsq;
        }
        return std::sqrt(best);
    }

    // ---------------------- Smoothed-from-signed helpers ----------------------

    // Given a signed distance value (positive inside, negative outside) produce a smooth result
    // using cubic smoothstep in a narrow band of half-width 'w'.
    // If w <= 0 returns the original signed value.
    static float SmoothSignedDistance_Smoothstep_FromSigned(float signedDist, float w)
    {
        if (w <= 0.0f) return signedDist;
        float sgn = (signedDist >= 0.0f) ? 1.0f : -1.0f;
        float a = std::fabs(signedDist);
        if (a >= w) return signedDist;
        float t = a / w; // in [0,1]
        float S = (3.0f * t * t - 2.0f * t * t * t); // cubic smoothstep
        float blended = S * w;
        return sgn * blended;
    }

    // Tanh-based smoothing that blends near origin controlled by alpha.
    // If alpha <= 0 returns original signedDist.
    // This version blends to exact distance outside small band to avoid saturation errors.
    static float SmoothSignedDistance_Tanh_FromSigned(float signedDist, float alpha)
    {
        if (alpha <= 0.0f) return signedDist;
        float a = std::fabs(signedDist);

        // If far away relative to alpha, return exact to avoid tanh saturation changing magnitude far from surface.
        if (a > 10.0f * alpha) return signedDist;

        // soft sign in (-1,1)
        float s = std::tanh(signedDist / alpha);

        // blend factor (0 near origin, ~1 far away)
        float blend = 0.5f + 0.5f * std::tanh((a - alpha) / (0.5f * alpha));

        // combine bounded tanh*alpha near 0 with exact signed distance far away
        return (1.0f - blend) * (s * alpha) + blend * signedDist;
    }



    // -------------------- BVH-accelerated nearest-triangle search --------------------

    // Squared distance from point p to AABB box
    static inline float AABBDistanceSq(const Vector3& p, const AABB& box)
    {
        float d2 = 0.0f;
        if (p.x < box.mn.x) { float t = box.mn.x - p.x; d2 += t*t; }
        else if (p.x > box.mx.x) { float t = p.x - box.mx.x; d2 += t*t; }

        if (p.y < box.mn.y) { float t = box.mn.y - p.y; d2 += t*t; }
        else if (p.y > box.mx.y) { float t = p.y - box.mx.y; d2 += t*t; }

        if (p.z < box.mn.z) { float t = box.mn.z - p.z; d2 += t*t; }
        else if (p.z > box.mx.z) { float t = p.z - box.mx.z; d2 += t*t; }

        return d2;
    }

    // BVH-accelerated closest-point search.
    // Returns true if a triangle was found (i.e., mesh non-empty).
    // Fills outClosest (closest point on mesh), outTri (triangle index), outDistSq (squared distance).
    bool ClosestPointBVH(const Vector3& p, Vector3& outClosest, int& outTri, float& outDistSq) const
    {
        // If no triangles, fail
        if (triangles.empty()) return false;

        // If BVH empty (not built), fall back to brute-force
        if (bvh.empty()) {
            // brute force search
            outDistSq = std::numeric_limits<float>::infinity();
            outTri = -1;
            Vector3 tmp;
            for (size_t i = 0; i < triangles.size(); ++i) {
                float d2 = PointTriangleDistanceSq(p, triangles[i], &tmp);
                if (d2 < outDistSq) {
                    outDistSq = d2;
                    outClosest = tmp;
                    outTri = static_cast<int>(i);
                }
            }
            return outTri != -1;
        }

        // Iterative stack traversal with pruning by AABB distance squared
        outDistSq = std::numeric_limits<float>::infinity();
        outTri = -1;

        // stack of node indices
        std::vector<int> stack;
        stack.reserve(64);
        stack.push_back(0); // root

        while (!stack.empty()) {
            int idx = stack.back();
            stack.pop_back();
            if (idx < 0) continue;
            const BVHNode& N = bvh[idx];

            // Compute lower-bound distance from p to this node's AABB:
            float boxDistSq = AABBDistanceSq(p, N.box);
            if (boxDistSq >= outDistSq) {
                // pruning: even best point inside this box can't beat current best
                continue;
            }

            if (N.count > 0) {
                // leaf: test all triangles
                for (int i = 0; i < N.count; ++i) {
                    int triIdx = primIndices[N.start + i];
                    Vector3 triClosest;
                    float d2 = PointTriangleDistanceSq(p, triangles[triIdx], &triClosest);
                    if (d2 < outDistSq) {
                        outDistSq = d2;
                        outClosest = triClosest;
                        outTri = triIdx;
                    }
                }
            } else {
                // internal: we push children in an order that prefers the child whose box is closer
                int left = N.left;
                int right = N.right;
                if (left < 0 && right < 0) continue;

                // If both children exist, compute their box distances and push the closer one first
                if (left >= 0 && right >= 0) {
                    float dl = AABBDistanceSq(p, bvh[left].box);
                    float dr = AABBDistanceSq(p, bvh[right].box);
                    // push farther first so closer one is processed sooner (LIFO stack)
                    if (dl < dr) {
                        if (dr < outDistSq) stack.push_back(right);
                        if (dl < outDistSq) stack.push_back(left);
                    } else {
                        if (dl < outDistSq) stack.push_back(left);
                        if (dr < outDistSq) stack.push_back(right);
                    }
                } else {
                    // one child only
                    if (left >= 0) {
                        if (AABBDistanceSq(p, bvh[left].box) < outDistSq) stack.push_back(left);
                    }
                    if (right >= 0) {
                        if (AABBDistanceSq(p, bvh[right].box) < outDistSq) stack.push_back(right);
                    }
                }
            }
        }

        return outTri != -1;
    }

    // Convenience wrapper: returns unsigned distance (sqrt) using BVH if available (fast).
    float UnsignedDistanceBVH(const Vector3& p) const
    {
        Vector3 closest;
        int tri;
        float d2;
        if (ClosestPointBVH(p, closest, tri, d2)) {
            return std::sqrt(d2);
        }
        // fallback
        return UnsignedDistanceBrute(p);
    }

    // ---------------------- SignedDistanceBVH with smoothing parameter ----------------------

    // Returns signed distance; positive = inside, negative = outside.
    // smoothingWidth:
    //   <= 0.0 : no smoothing, exact signed distance (uses BVH nearest if available).
    //   > 0.0  : smoothing parameter. If useTanh==false -> cubic smoothstep narrow-band using smoothingWidth.
    //            If useTanh==true  -> tanh-based smoothing with alpha = smoothingWidth.
    float SignedDistanceBVH(const Vector3& p, float smoothingWidth = 0.0f, bool useTanh = false) const
    {
        // Compute unsigned distance using BVH if available (fast).
        float ud = 0.0f;
        if (!bvh.empty()) ud = UnsignedDistanceBVH(p);
        else ud = UnsignedDistanceBrute(p);

        // Determine sign (inside => positive)
        bool inside = IsPointInsideObject(p);
        float signedDist = inside ? ud : -ud;

        // Apply smoothing if requested
        if (smoothingWidth <= 0.0f) return signedDist;

        if (useTanh) {
            return SmoothSignedDistance_Tanh_FromSigned(signedDist, smoothingWidth);
        } else {
            return SmoothSignedDistance_Smoothstep_FromSigned(signedDist, smoothingWidth);
        }
    }


    // ----------------- end helpers --------------------------------

    // Constructor
    STLObject() { ResetTransform(); }
};

}  // namespace StlGeom

#endif // STL_READER_H_
