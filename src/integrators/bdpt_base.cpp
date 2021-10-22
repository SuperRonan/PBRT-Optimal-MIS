
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
        path[0].pdfFwd = pdfPos;
        // With float, it might be zero.
        // I did not observe this problem with double.
        if (pdfPos == 0)
            return 0;
        if (pdfDir == 0)
            return 1;
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
        Vertex& lv = path[0];
        if (lv.IsInfiniteLight()) {
            // Set spatial density of _path[1]_ for infinite area light
            if (lv.ei.light->flags & (int)LightFlags::DeltaDirection) // Directional light
            {
                lv.pdfFwd = pdfDir;
                lv.pdfRev = 0;
            }
            else // Environement map
            {
                lv.pdfFwd = InfiniteLightDensity(scene, lightDistr, lightToIndex, ray.d);
            }
            if (nVertices > 0) {
                path[1].pdfFwd = pdfPos;
                if (path[1].IsOnSurface())
                    path[1].pdfFwd *= AbsDot(ray.d, path[1].ng());
            }
        }
        return nVertices + 1;
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

    Float MISWeight(const Scene& scene, Vertex* lightVertices,
        Vertex* cameraVertices, Vertex& sampled, int s, int t,
        const Distribution1D& lightPdf,
        const std::unordered_map<const Light*, size_t>& lightToIndex, Float s1Pdf) {
        if (s + t == 2) return 1;
        Float sumRi = 0;

        // Temporarily update vertex properties for current strategy

        // Look up connection vertices and their predecessors
        Vertex* qs = s > 0 ? &lightVertices[s - 1] : nullptr,
            * pt = t > 0 ? &cameraVertices[t - 1] : nullptr,
            * qsMinus = s > 1 ? &lightVertices[s - 2] : nullptr,
            * ptMinus = t > 1 ? &cameraVertices[t - 2] : nullptr;

        // Update sampled vertex for $s=1$ or $t=1$ strategy
        ScopedAssignment<Vertex> a1;
        if (s == 1)
            a1 = { qs, sampled };
        else if (t == 1)
            a1 = { pt, sampled };

        // Mark connection vertices as non-degenerate
        ScopedAssignment<bool> a2, a3;
        if (pt) a2 = { &pt->delta, false };
        if (qs) a3 = { &qs->delta, false };

        // Update reverse density of vertex $\pt{}_{t-1}$
        ScopedAssignment<Float> a4;
        if (pt)
            a4 = { &pt->pdfRev, s > 0 ? qs->Pdf(scene, qsMinus, *pt)
                                     : pt->PdfLightOrigin(scene, *ptMinus, lightPdf,
                                                          lightToIndex) };

        // Update reverse density of vertex $\pt{}_{t-2}$
        ScopedAssignment<Float> a5;
        if (ptMinus)
            a5 = { &ptMinus->pdfRev, s > 0 ? pt->Pdf(scene, qs, *ptMinus)
                                          : pt->PdfLight(scene, *ptMinus) };

        // Update reverse density of vertices $\pq{}_{s-1}$ and $\pq{}_{s-2}$
        ScopedAssignment<Float> a6;
        if (qs) a6 = { &qs->pdfRev, pt->Pdf(scene, ptMinus, *qs) };
        ScopedAssignment<Float> a7;
        if (qsMinus) a7 = { &qsMinus->pdfRev, qs->Pdf(scene, pt, *qsMinus) };

        // Consider hypothetical connection strategies along the camera subpath
        Float ri = 1;
        for (int i = t - 1; i > 0; --i) {
            ri *=
                cameraVertices[i].pdfRev / cameraVertices[i].pdfFwd;
            Float actualRi = (s == 0 && i == t - 1) ?
                ri * s1Pdf / cameraVertices[i].pdfRev : ri;
            if (!cameraVertices[i].delta && !cameraVertices[i - 1].delta)
                sumRi += actualRi;
        }

        // Consider hypothetical connection strategies along the light subpath
        ri = 1;
        for (int i = s - 1; i >= 0; --i) {
            ri *= lightVertices[i].pdfRev / lightVertices[i].pdfFwd;
            bool deltaLightvertex = i > 0 ? lightVertices[i - 1].delta
                : lightVertices[0].IsDeltaLight();
            Float actualRi = (i == 1) ?
                ri * s1Pdf / lightVertices[0].pdfFwd : ri;
            if (!lightVertices[i].delta && !deltaLightvertex) sumRi += actualRi;
        }
        Float actualRi = (s == 1) ? s1Pdf / lightVertices[0].pdfFwd : 1.0;
        return actualRi / (actualRi + sumRi);
    }


    Spectrum ConnectBDPT(
        const Scene& scene, Vertex* lightVertices, Vertex* cameraVertices, int s,
        int t, const Distribution1D& lightDistr,
        const std::unordered_map<const Light*, size_t>& lightToIndex,
        const Camera& camera, Sampler& sampler, Point2f* pRaster,
        Float* misWeightPtr) {
        ProfilePhase _(Prof::BDPTConnectSubpaths);
        Spectrum L(0.f);
        // Ignore invalid connections related to infinite area lights
        if (t > 1 && s != 0 && cameraVertices[t - 1].type == VertexType::Light)
            return Spectrum(0.f);

        Float s1Pdf;
        // Perform connection and write contribution to _L_
        Vertex sampled;
        if (s == 0) {
            // Interpret the camera subpath as a complete path
            const Vertex& pt = cameraVertices[t - 1];
            if (pt.IsLight()) {
                L = pt.Le(scene, cameraVertices[t - 2]) * pt.beta;
                s1Pdf = cameraVertices[t - 1].PdfResampledLight(scene, cameraVertices[t - 2], lightDistr, lightToIndex);
            }
            DCHECK(!L.HasNaNs());
        }
        else if (t == 1) {
            // Sample a point on the camera and connect it to the light subpath
            const Vertex& qs = lightVertices[s - 1];
            if (qs.IsConnectible()) {
                VisibilityTester vis;
                Vector3f wi;
                Float pdf;
                Spectrum Wi = camera.Sample_Wi(qs.GetInteraction(), sampler.Get2D(),
                    &wi, &pdf, pRaster, &vis);
                if (pdf > 0 && !Wi.IsBlack()) {
                    // Initialize dynamically sampled vertex and _L_ for $t=1$ case
                    sampled = Vertex::CreateCamera(&camera, vis.P1(), Wi / pdf);
                    sampled.pdfFwd = pdf;
                    L = qs.beta * qs.f(sampled, TransportMode::Importance) * sampled.beta;
                    if (qs.IsOnSurface()) L *= AbsDot(wi, qs.ns());
                    DCHECK(!L.HasNaNs());
                    // Only check visibility after we know that the path would
                    // make a non-zero contribution.
                    if (!L.IsBlack()) L *= vis.Tr(scene, sampler);
                }
            }
        }
        else if (s == 1) {
            // Sample a point on a light and connect it to the camera subpath
            const Vertex& pt = cameraVertices[t - 1];
            if (pt.IsConnectible()) {
                Float lightPdf;
                VisibilityTester vis;
                Vector3f wi;
                Float pdfSolidAngle;
                int lightNum =
                    lightDistr.SampleDiscrete(sampler.Get1D(), &lightPdf);
                const std::shared_ptr<Light>& light = scene.lights[lightNum];
                Spectrum lightWeight = light->Sample_Li(
                    pt.GetInteraction(), sampler.Get2D(), &wi, &pdfSolidAngle, &vis);
                if (pdfSolidAngle > 0 && !lightWeight.IsBlack()) {
                    EndpointInteraction ei(vis.P1(), light.get());
                    sampled =
                        Vertex::CreateLight(ei, lightWeight / (pdfSolidAngle * lightPdf), 0);
                    sampled.pdfFwd =
                        sampled.PdfLightOrigin(scene, pt, lightDistr, lightToIndex);
                    Float pdfArea = pt.ConvertDensity(pdfSolidAngle, sampled);
                    s1Pdf = pdfArea * lightPdf;
                    L = pt.beta * pt.f(sampled, TransportMode::Radiance) * sampled.beta;
                    if (pt.IsOnSurface()) L *= AbsDot(wi, pt.ns());
                    // Only check visibility if the path would carry radiance.
                    if (!L.IsBlack()) L *= vis.Tr(scene, sampler);
                }
            }
        }
        else {
            // Handle all other bidirectional connection cases
            const Vertex& qs = lightVertices[s - 1], & pt = cameraVertices[t - 1];
            if (qs.IsConnectible() && pt.IsConnectible()) {
                L = qs.beta * qs.f(pt, TransportMode::Importance) * pt.f(qs, TransportMode::Radiance) * pt.beta;
                VLOG(2) << "General connect s: " << s << ", t: " << t <<
                    " qs: " << qs << ", pt: " << pt << ", qs.f(pt): " << qs.f(pt, TransportMode::Importance) <<
                    ", pt.f(qs): " << pt.f(qs, TransportMode::Radiance) << ", G: " << G(scene, sampler, qs, pt) <<
                    ", dist^2: " << DistanceSquared(qs.p(), pt.p());
                if (!L.IsBlack()) L *= G(scene, sampler, qs, pt);
            }
        }
        if (s >= 2)
            s1Pdf = lightVertices[0].PdfResampledLight(scene, lightVertices[1], lightDistr, lightToIndex);

        // Compute MIS weight for connection strategy
        Float misWeight =
            L.IsBlack() ? 0.f : MISWeight(scene, lightVertices, cameraVertices,
                sampled, s, t, lightDistr, lightToIndex, s1Pdf);
        VLOG(2) << "MIS weight for (s,t) = (" << s << ", " << t << ") connection: "
            << misWeight;
        DCHECK(!std::isnan(misWeight));
        L *= misWeight;
        if (misWeightPtr) *misWeightPtr = misWeight;
        return L;
    }
}