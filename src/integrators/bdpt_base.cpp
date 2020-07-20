
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

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

#include "bdpt_base.h"

namespace pbrt
{
    Float CorrectShadingNormal(const SurfaceInteraction& isect, const Vector3f& wo,
        const Vector3f& wi, TransportMode mode) {
        if (mode == TransportMode::Importance) {
            Float num = AbsDot(wo, isect.shading.n) * AbsDot(wi, isect.n);
            Float denom = AbsDot(wo, isect.n) * AbsDot(wi, isect.shading.n);
            // wi is occasionally perpendicular to isect.shading.n; this is
            // fine, but we don't want to return an infinite or NaN value in
            // that case.
            if (denom == 0) return 0;
            return num / denom;
        }
        else
            return 1;
    }

    int GenerateCameraSubpath(const Scene& scene, Sampler& sampler,
        MemoryArena& arena, int maxDepth,
        const Camera& camera, const Point2f& pFilm,
        Vertex* path) {
        if (maxDepth == 0) return 0;
        ProfilePhase _(Prof::BDPTGenerateSubpath);
        // Sample initial ray for camera subpath
        CameraSample cameraSample;
        cameraSample.pFilm = pFilm;
        cameraSample.time = sampler.Get1D();
        cameraSample.pLens = sampler.Get2D();
        RayDifferential ray;
        Spectrum beta = camera.GenerateRayDifferential(cameraSample, &ray);
        ray.ScaleDifferentials(1 / std::sqrt(sampler.samplesPerPixel));

        // Generate first vertex on camera subpath and start random walk
        Float pdfPos, pdfDir;
        path[0] = Vertex::CreateCamera(&camera, ray, beta);
        camera.Pdf_We(ray, &pdfPos, &pdfDir);
        VLOG(2) << "Starting camera subpath. Ray: " << ray << ", beta " << beta
            << ", pdfPos " << pdfPos << ", pdfDir " << pdfDir;
        return RandomWalk(scene, ray, sampler, arena, beta, pdfDir, maxDepth - 1,
            TransportMode::Radiance, path + 1) +
            1;
    }

    int GenerateLightSubpath(
        const Scene& scene, Sampler& sampler, MemoryArena& arena, int maxDepth,
        Float time, const Distribution1D& lightDistr,
        const std::unordered_map<const Light*, size_t>& lightToIndex,
        Vertex* path) {
        if (maxDepth == 0) return 0;
        ProfilePhase _(Prof::BDPTGenerateSubpath);
        // Sample initial ray for light subpath
        Float lightPdf;
        int lightNum = lightDistr.SampleDiscrete(sampler.Get1D(), &lightPdf);
        const std::shared_ptr<Light>& light = scene.lights[lightNum];
        RayDifferential ray;
        Normal3f nLight;
        Float pdfPos, pdfDir;
        Spectrum Le = light->Sample_Le(sampler.Get2D(), sampler.Get2D(), time, &ray,
            &nLight, &pdfPos, &pdfDir);
        if (pdfPos == 0 || pdfDir == 0 || Le.IsBlack()) return 0;

        // Generate first vertex on light subpath and start random walk
        path[0] =
            Vertex::CreateLight(light.get(), ray, nLight, Le, pdfPos * lightPdf);
        Spectrum beta = Le * AbsDot(nLight, ray.d) / (lightPdf * pdfPos * pdfDir);
        VLOG(2) << "Starting light subpath. Ray: " << ray << ", Le " << Le <<
            ", beta " << beta << ", pdfPos " << pdfPos << ", pdfDir " << pdfDir;
        int nVertices =
            RandomWalk(scene, ray, sampler, arena, beta, pdfDir, maxDepth - 1,
                TransportMode::Importance, path + 1);

        // Correct subpath sampling densities for infinite area lights
        if (path[0].IsInfiniteLight()) {
            // Set spatial density of _path[1]_ for infinite area light
            if (nVertices > 0) {
                path[1].pdfFwd = pdfPos;
                if (path[1].IsOnSurface())
                    path[1].pdfFwd *= AbsDot(ray.d, path[1].ng());
            }

            // Set spatial density of _path[0]_ for infinite area light
            path[0].pdfFwd =
                InfiniteLightDensity(scene, lightDistr, lightToIndex, ray.d);
        }
        return nVertices + 1;
    }

    int RandomWalk(const Scene& scene, RayDifferential ray, Sampler& sampler,
        MemoryArena& arena, Spectrum beta, Float pdf, int maxDepth,
        TransportMode mode, Vertex* path) {
        if (maxDepth == 0) return 0;
        int bounces = 0;
        // Declare variables for forward and reverse probability densities
        Float pdfFwd = pdf, pdfRev = 0;
        while (true) {
            // Attempt to create the next subpath vertex in _path_
            MediumInteraction mi;

            VLOG(2) << "Random walk. Bounces " << bounces << ", beta " << beta <<
                ", pdfFwd " << pdfFwd << ", pdfRev " << pdfRev;
            // Trace a ray and sample the medium, if any
            SurfaceInteraction isect;
            bool foundIntersection = scene.Intersect(ray, &isect);
            if (ray.medium) beta *= ray.medium->Sample(ray, sampler, arena, &mi);
            if (beta.IsBlack()) break;
            Vertex& vertex = path[bounces], & prev = path[bounces - 1];
            if (mi.IsValid()) {
                // Record medium interaction in _path_ and compute forward density
                vertex = Vertex::CreateMedium(mi, beta, pdfFwd, prev);
                if (++bounces >= maxDepth) break;

                // Sample direction and compute reverse density at preceding vertex
                Vector3f wi;
                pdfFwd = pdfRev = mi.phase->Sample_p(-ray.d, &wi, sampler.Get2D());
                ray = mi.SpawnRay(wi);
            }
            else {
                // Handle surface interaction for path generation
                if (!foundIntersection) {
                    // Capture escaped rays when tracing from the camera
                    if (mode == TransportMode::Radiance) {
                        vertex = Vertex::CreateLight(EndpointInteraction(ray), beta,
                            pdfFwd);
                        ++bounces;
                    }
                    break;
                }

                // Compute scattering functions for _mode_ and skip over medium
                // boundaries
                isect.ComputeScatteringFunctions(ray, arena, true, mode);
                if (!isect.bsdf) {
                    ray = isect.SpawnRay(ray.d);
                    continue;
                }

                // Initialize _vertex_ with surface intersection information
                vertex = Vertex::CreateSurface(isect, beta, pdfFwd, prev);
                if (++bounces >= maxDepth) break;

                // Sample BSDF at current vertex and compute reverse probability
                Vector3f wi, wo = isect.wo;
                BxDFType type;
                Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdfFwd,
                    BSDF_ALL, &type);
                VLOG(2) << "Random walk sampled dir " << wi << " f: " << f <<
                    ", pdfFwd: " << pdfFwd;
                if (f.IsBlack() || pdfFwd == 0.f) break;
                beta *= f * AbsDot(wi, isect.shading.n) / pdfFwd;
                VLOG(2) << "Random walk beta now " << beta;
                if (type & BSDF_SPECULAR) {
                    vertex.delta = true;
                    pdfRev = pdfFwd; // This heavily assumes that the specular pdf is symmetric
                }
                else
                    pdfRev = isect.bsdf->Pdf(wi, wo, BSDF_ALL);
                beta *= CorrectShadingNormal(isect, wo, wi, mode);
                VLOG(2) << "Random walk beta after shading normal correction " << beta;
                ray = isect.SpawnRay(wi);
            }

            // Compute reverse area density at preceding vertex
            prev.pdfRev = vertex.ConvertDensity(pdfRev, prev);
        }
        return bounces;
    }

    Spectrum G(const Scene& scene, Sampler& sampler, const Vertex& v0,
        const Vertex& v1) {
        Vector3f d = v0.p() - v1.p();
        Float g = 1 / d.LengthSquared();
        d *= std::sqrt(g);
        if (v0.IsOnSurface()) g *= AbsDot(v0.ns(), d);
        if (v1.IsOnSurface()) g *= AbsDot(v1.ns(), d);
        VisibilityTester vis(v0.GetInteraction(), v1.GetInteraction());
        return g * vis.Tr(scene, sampler);
    }
}