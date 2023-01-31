#pragma once

#include <core/shape.h>

#include <vector>
#include <memory>

namespace lightfold {

    struct TriangleMesh {
        // TriangleMesh Public Methods
        TriangleMesh(const Transform& ObjectToWorld, int nTriangles, const int* vertexIndices,
            int nVertices, const Point3f* P, const Tangent3f* S, const Normal3f* N, const Point2f* UV)
            : nTriangles(nTriangles), nVertices(nVertices),
            vertexIndices(vertexIndices, vertexIndices + 3 * nTriangles) {
            p.reset(new Point3f[nVertices]);
            for (int i = 0; i < nVertices; ++i)
                p[i] = ObjectToWorld(P[i]);
            if (UV) {
                uv.reset(new Point2f[nVertices]);
                memcpy(uv.get(), UV, nVertices * sizeof(Point2f));
            }
            if (N) {
                n.reset(new Normal3f[nVertices]);
                for (int i = 0; i < nVertices; ++i)
                    n[i] = ObjectToWorld(N[i]);
            }
            if (S) {
                s.reset(new Tangent3f[nVertices]);
                for (int i = 0; i < nVertices; ++i)
                    s[i] = ObjectToWorld(S[i]);
            }
        }

        // TriangleMesh Public Data
        const int nTriangles, nVertices;
        std::vector<int> vertexIndices;
        std::unique_ptr<Point3f[]> p;
        std::unique_ptr<Normal3f[]> n;
        std::unique_ptr<Tangent3f[]> s;
        std::unique_ptr<Point2f[]> uv;
    };

    std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
        const Transform* ObjectToWorld, const Transform* WorldToObject,
        bool reverseOrientation, int nTriangles, const int* vertexIndices,
        int nVertices, const Point3f* p, const Tangent3f* s, const Normal3f* n,
        const Point2f* uv);

    class Triangle : public Shape {
    public:
        // Triangle Public Methods
        Triangle(const Transform* ObjectToWorld, const Transform* WorldToObject,
            bool reverseOrientation,
            const std::shared_ptr<TriangleMesh>& mesh, int triNumber)
            : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
            mesh(mesh) {
            v = &mesh->vertexIndices[3 * triNumber];
        }

        Bounds3f ObjectBound() const;
        Bounds3f WorldBound() const;

        bool Intersect(const Ray& ray, float* tHit, SurfaceInteraction* isect) const;
        //bool IntersectP(const Ray& ray) const;

    private:
        // Triangle Private Methods
        void GetUVs(Point2f uv[3]) const {
            if (mesh->uv) {
                uv[0] = mesh->uv[v[0]];
                uv[1] = mesh->uv[v[1]];
                uv[2] = mesh->uv[v[2]];
            }
            else {
                uv[0] = Point2f(0, 0);
                uv[1] = Point2f(1, 0);
                uv[2] = Point2f(1, 1);
            }
        }

        // Triangle Private Data
        std::shared_ptr<TriangleMesh> mesh;
        const int* v;
    };

} // namespace lightfold