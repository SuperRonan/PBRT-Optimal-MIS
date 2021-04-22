
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

 // integrators/bdpt.cpp*
#include "integrators/bdpt_opti.h"
#include "film.h"
#include "filters/box.h"
#include "integrator.h"
#include "lightdistrib.h"
#include "paramset.h"
#include "progressreporter.h"
#include "sampler.h"
#include "stats.h"


// Uncorrelated BDPT
// Draw a new camera and light sub path for every connection strategy
// Very ineficient, for research purpose only
//#define UBDPT

namespace pbrt {

    STAT_PERCENT("Integrator/Zero-radiance paths", zeroRadiancePaths, totalPaths);
    STAT_INT_DISTRIBUTION("Integrator/Path length", pathLength);


    // OBDPT Utility Functions
    // Fills the weights vector with the balance weights
    void BalanceWeights(const Scene& scene, Vertex* lightVertices,
        Vertex* cameraVertices, Vertex& sampled, const int s, const int t,
        const Distribution1D& lightPdf,
        const std::unordered_map<const Light*, size_t>& lightToIndex, Float s1Pdf, 
        Float * weights, bool & sparseZero) {
        DCHECK_NE(s + t, 2);
        Float sumRi = 0;
        Float * const& ratios = weights;

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
        if (pt) {
            Float pdf = ((s == 0) ?
                pt->PdfLightOrigin(scene, *ptMinus, lightPdf, lightToIndex) :
                qs->Pdf(scene, qsMinus, *pt));
            if (pdf == 0 && s == 0)
            {
                // pt is on the env map, but the env map is dark, so not really a light source
                sparseZero = true;
                return;
            }
            a4 = { &pt->pdfRev, pdf};
        }

        // Update reverse density of vertex $\pt{}_{t-2}$
        ScopedAssignment<Float> a5;
        if (ptMinus && ptMinus->type != VertexType::Camera)
        {
            Float pdf = (s == 0) ?
                pt->PdfLight(scene, *ptMinus) :
                pt->Pdf(scene, qs, *ptMinus);
            a5 = { &ptMinus->pdfRev, pdf };
        }
        
        // Update reverse density of vertices $\pq{}_{s-1}$ and $\pq{}_{s-2}$
        ScopedAssignment<Float> a6;
        if (qs) a6 = { &qs->pdfRev, pt->Pdf(scene, ptMinus, *qs) };
        ScopedAssignment<Float> a7;
        if (qsMinus) a7 = { &qsMinus->pdfRev, qs->Pdf(scene, pt, *qsMinus) };

        // Consider hypothetical connection strategies along the camera subpath
        Float ri = 1;
        for (int i = t - 1; i > 0; --i) {
            ri *= cameraVertices[i].pdfRev / cameraVertices[i].pdfFwd;
            Float actualRi = (s == 0 && i == t - 1) ?
                ri * s1Pdf / cameraVertices[i].pdfRev : ri;
            if (cameraVertices[i].delta | cameraVertices[i - 1].delta)
                actualRi = 0;
            ratios[(s+t)-i] = actualRi;
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
            if (lightVertices[i].delta | deltaLightvertex)
                actualRi = 0;
            sumRi += actualRi;
            ratios[i] = actualRi;
        }
        Float actualRi = (s == 1) ? s1Pdf / lightVertices[0].pdfFwd : 1.0;
        ratios[s] = actualRi;
        sumRi += actualRi;

        DCHECK_GT(sumRi, 0);

        for (int i = 0; i < (s+t); ++i)
            weights[i] = ratios[i] / sumRi ;

    }

    // BDPT Method Definitions
    inline int BufferIndex(int s, int t) {
        int above = s + t - 2;
        return s + above * (5 + above) / 2;
    }

    void OBDPTIntegrator::Render(const Scene& scene) {
        std::unique_ptr<LightDistribution> lightDistribution =
            CreateLightSampleDistribution(lightSampleStrategy, scene);

        std::cout << "Rendering with setting: \n";
        std::cout << "Conservative: " << conservative << std::endl;
#ifdef UBDPT
        std::cout << "Uncorellated: 1" << std::endl;
#else   
        std::cout << "Uncorellated: 0" << std::endl;
#endif

        // Compute a reverse mapping from light pointers to offsets into the
        // scene lights vector (and, equivalently, offsets into
        // lightDistr). Added after book text was finalized; this is critical
        // to reasonable performance with 100s+ of light sources.
        std::unordered_map<const Light*, size_t> lightToIndex;
        for (size_t i = 0; i < scene.lights.size(); ++i)
            lightToIndex[scene.lights[i].get()] = i;

        // Partition the image into tiles
        Film* film = camera->film;
        const Bounds2i sampleBounds = film->GetSampleBounds();
        const Vector2i sampleExtent = sampleBounds.Diagonal();
        const int tileSize = 16;
        const int nXTiles = (sampleExtent.x + tileSize - 1) / tileSize;
        const int nYTiles = (sampleExtent.y + tileSize - 1) / tileSize;

        const size_t image_seed = size_t(nXTiles * nYTiles) * size_t(sampler->samplesPerPixel) * this->seed;

        // Render and write the output image to disk
        if (scene.lights.size() > 0) {
            // Pre allocate estimators
            const bool USE_ROW_MAJOR = true;
            std::vector<MIS::ImageEstimator<Spectrum, Float, USE_ROW_MAJOR>*> estimators;
            for (int depth = std::max(1, minDepth); depth <= maxOptiDepth; ++depth)
            {
                int N = depth + 2;
                int width = sampleExtent.x;
                int height = sampleExtent.y;
                estimators.emplace_back(MIS::createImageEstimator<Spectrum, Float, USE_ROW_MAJOR>(heuristic, N, width, height));
                estimators.back()->setSampleForTechnique(N - 1, width * height);
            }
            {
                ProgressReporter reporter(nXTiles * nYTiles, "Drawing Samples");

                ParallelFor2D([&](const Point2i tile) {
                    // Render a single tile using BDPT
                    MemoryArena arena;
                    int seed = tile.y * nXTiles + tile.x + image_seed;
                    std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);
                    int x0 = sampleBounds.pMin.x + tile.x * tileSize;
                    int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
                    int y0 = sampleBounds.pMin.y + tile.y * tileSize;
                    int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
                    Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
                    LOG(INFO) << "Starting image tile " << tileBounds;

                    std::unique_ptr<FilmTile> filmTile =
                        camera->film->GetFilmTile(tileBounds);
                    for (Point2i pPixel : tileBounds) {
                        tileSampler->StartPixel(pPixel);
                        if (!InsideExclusive(pPixel, pixelBounds))
                            continue;
                        do {
                            // Generate a single sample using BDPT
                            Point2f pFilm = (Point2f)pPixel + tileSampler->Get2D();

                            // Trace the camera subpath
                            Vertex* cameraVertices = arena.Alloc<Vertex>(maxDepth + 2);
                            Vertex* lightVertices = arena.Alloc<Vertex>(maxDepth + 1);
                            Float* balanceWeights = arena.Alloc<Float>(maxDepth + 2);

                            // Get a distribution for sampling the light at the
                            // start of the light subpath. Because the light path
                            // follows multiple bounces, basing the sampling
                            // distribution on any of the vertices of the camera
                            // path is unlikely to be a good strategy. We use the
                            // PowerLightDistribution by default here, which
                            // doesn't use the point passed to it.
                            const Distribution1D* lightDistr =
                                lightDistribution->Lookup(cameraVertices[0].p());

#ifndef UBDPT
                            int nCamera = GenerateCameraSubpath(
                                scene, *tileSampler, arena, maxDepth + 2, *camera,
                                pFilm, cameraVertices);
                            // Now trace the light subpath
                            int nLight = GenerateLightSubpath(
                                scene, *tileSampler, arena, maxDepth + 1,
                                cameraVertices[0].time(), *lightDistr, lightToIndex,
                                lightVertices);
#else
                            int nCamera = maxDepth + 2;
                            int nLight = maxDepth + 1;
#endif

                            // Execute all BDPT connection strategies
                            Spectrum estimate(0.f);
                            Spectrum L = 0;
                            for (int t = 1; t <= nCamera; ++t) {
                                for (int s = 0; s <= nLight; ++s) {
                                    int depth = t + s - 2;
                                    if ((s == 1 && t == 1) || depth < minDepth ||
                                        depth > maxDepth)
                                        continue;
                                    bool fallBackBalance = (depth == 0) || (depth > maxOptiDepth);
#ifdef UBDPT
                                    int nCamera_generated = GenerateCameraSubpath(
                                        scene, *tileSampler, arena, t, *camera,
                                        pFilm, cameraVertices);
                                    if (nCamera_generated != t)  continue;
                                    // Now trace the light subpath
                                    int nLight_generated = GenerateLightSubpath(
                                        scene, *tileSampler, arena, s,
                                        cameraVertices[0].time(), *lightDistr, lightToIndex,
                                        lightVertices);
                                    if (nLight_generated != s)  continue;
#endif

                                    Point2f pFilmNew = pFilm;
                                    if (fallBackBalance)
                                    {
                                        estimate = ConnectBDPT(scene, lightVertices, cameraVertices, s, t,
                                            *lightDistr, lightToIndex, *camera, *tileSampler,
                                            &pFilmNew);
                                        ++totalPaths;
                                        if (estimate.IsBlack()) ++zeroRadiancePaths;
                                        ReportValue(pathLength, s + t - 2);
                                        if (!estimate.IsBlack())
                                            if (t != 1)
                                                L += estimate;
                                            else
                                                film->AddSplat(pFilmNew, estimate / sampler->samplesPerPixel);
                                    }
                                    else
                                    {
                                        bool sparseZero = true;
                                        if(conservative)
                                            estimate = ConnectOBDPT<true>(
                                                scene, lightVertices, cameraVertices, s, t,
                                                *lightDistr, lightToIndex, *camera, *tileSampler,
                                                &pFilmNew, balanceWeights, sparseZero);
                                        else
                                            estimate = ConnectOBDPT<false>(
                                                scene, lightVertices, cameraVertices, s, t,
                                                *lightDistr, lightToIndex, *camera, *tileSampler,
                                                &pFilmNew, balanceWeights, sparseZero);

                                        if (!sparseZero)
                                        {
                                            // This is maybe not very clever, since the inverse is done in the addEstimate...
                                            Float u = pFilmNew.x / sampleExtent.x;
                                            Float v = pFilmNew.y / sampleExtent.y;
                                            int index = depth - std::max(minDepth, 1);
                                            auto* estimator = estimators[index];
                                            estimator->addEstimate(estimate, balanceWeights, s, u, v);
                                        }
                                    }
                                } // for s
                            } // for s
                            //VLOG(2) << "Add film sample pFilm: " << pFilm << ", L: " << L <<
                            //    ", (y: " << L.y() << ")";
                            filmTile->AddSample(pFilm, L);
                            arena.Reset();
                        } while (tileSampler->StartNextSample());
                    }
                    film->MergeFilmTile(std::move(filmTile));
                    reporter.Update();
                    LOG(INFO) << "Finished image tile " << tileBounds;
                    }, Point2i(nXTiles, nYTiles));
                reporter.Done();
            }
            if(!estimators.empty())
            {
                // Make sure that OpenMP is enabled to launch the solving in parallel
                // It should be enabled in the cmakefiles 
                MIS::Parallel::init();
                MIS::Parallel::setNumThreads(MaxThreadIndex());
                ProgressReporter reporter(estimators.size(), "Solving");
                std::vector<Spectrum> buffer(sampleExtent.x * sampleExtent.y, Spectrum(0));
                for (auto * estimator : estimators)
                {
                    estimator->solve(buffer.data(), sampler->samplesPerPixel);
                    reporter.Update();
                }
                reporter.Done();

                estimators[0]->loopThroughImage([&film, &buffer, &estimators](int i, int j)
                    {
                        Point2f p(i + 0.5, j + 0.5);
                        film->AddSplat(p, buffer[estimators[0]->pixelTo1D(i, j)]);
                    });
            }
            for (int i = 0; i < estimators.size(); ++i)
                delete estimators[i], estimators[i] = nullptr;

        } // if scene has lights
        

        film->WriteImage();
    }

    // This function also needs to prove that the sample is not sparse zero
    template <bool CONSERVATIVE>
    Spectrum ConnectOBDPT(
        const Scene& scene, Vertex* lightVertices, Vertex* cameraVertices, int s,
        int t, const Distribution1D& lightDistr,
        const std::unordered_map<const Light*, size_t>& lightToIndex,
        const Camera& camera, Sampler& sampler, Point2f* pRaster,
        Float* balance_weights, bool& sparseZero) 
    {
        ProfilePhase _(Prof::BDPTConnectSubpaths);
        Spectrum L(0.f);
        // Ignore invalid connections related to infinite area lights
        if (t > 1 && s != 0 && cameraVertices[t - 1].type == VertexType::Light)
        {
            return 0;
        }
        DCHECK_NE(s + t, 2);
        Float s1Pdf;
        // Perform connection and write contribution to _L_
        Vertex sampled;
        if (s == 0) {
            // Interpret the camera subpath as a complete path
            const Vertex& pt = cameraVertices[t - 1];
            if (pt.IsLight()) {
                L = pt.Le(scene, cameraVertices[t - 2]) * pt.beta;
                s1Pdf = cameraVertices[t - 1].PdfResampledLight(scene, cameraVertices[t - 2], lightDistr, lightToIndex);
                sparseZero = false;
            }
            else
                return 0;
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
                    // make a non-zero contribution. (for non conservative estimation)
                    bool check_vis = CONSERVATIVE || !L.IsBlack();
                    if (check_vis)
                    {
                        Spectrum V = vis.Tr(scene, sampler);
                        L *= V;
                        sparseZero = V.IsBlack();
                    }
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
                    // Only check visibility if the path would carry radiance. (for conservative estimation)
                    bool check_vis = CONSERVATIVE || !L.IsBlack();
                    if (check_vis)
                    {
                        Spectrum V = vis.Tr(scene, sampler);
                        L *= V;
                        sparseZero = V.IsBlack();
                    }
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
                bool check_vis = CONSERVATIVE || !L.IsBlack();
                if (check_vis)
                {
                    Spectrum geometry = G(scene, sampler, qs, pt);
                    L *= geometry;
                    sparseZero = geometry.IsBlack();
                }
            }
        }
        if (s >= 2)
            s1Pdf = lightVertices[0].PdfResampledLight(scene, lightVertices[1], lightDistr, lightToIndex);

        ++totalPaths;
        if (L.IsBlack()) ++zeroRadiancePaths;
        ReportValue(pathLength, s + t - 2);


        if (!sparseZero)
        {
            BalanceWeights(scene, lightVertices, cameraVertices,
                sampled, s, t, lightDistr, lightToIndex, s1Pdf, balance_weights, sparseZero);
        }
        if (sparseZero) L = 0;
        return L;
    }

    OBDPTIntegrator* CreateOBDPTIntegrator(const ParamSet& params,
        std::shared_ptr<Sampler> sampler,
        std::shared_ptr<const Camera> camera) {
        int maxDepth = params.FindOneInt("maxdepth", 5);
        int minDepth = params.FindOneInt("mindepth", 0);
        int maxOptiDepth = params.FindOneInt("maxoptidepth", maxDepth);

        
        int np;
        const int* pb = params.FindInt("pixelbounds", &np);
        Bounds2i pixelBounds = camera->film->GetSampleBounds();
        if (pb) {
            if (np != 4)
                Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                    np);
            else {
                pixelBounds = Intersect(pixelBounds,
                    Bounds2i{ {pb[0], pb[2]}, {pb[1], pb[3]} });
                if (pixelBounds.Area() == 0)
                    Error("Degenerate \"pixelbounds\" specified.");
            }
        }

        std::string h_name = params.FindOneString("heuristic", "balance");
        MIS::Heuristic h;
        if (h_name == "balance")
            h = MIS::Heuristic::Balance;
        else if (h_name == "power")
            h = MIS::Heuristic::Power;
        else if (h_name == "naive")
            h = MIS::Heuristic::Naive;
        else if (h_name == "cutoff")
            h = MIS::Heuristic::CutOff;
        else if (h_name == "maximum")
            h = MIS::Heuristic::Maximum;
        else if (h_name == "direct")
            h = MIS::Heuristic::Direct;
        else
            Error(std::string(std::string("Heuristic ") + h_name + " is not recognized!").c_str());

        bool conservative = params.FindOneBool("conservative", false);

        std::string lightStrategy = params.FindOneString("lightsamplestrategy", "power");

        size_t seed = params.FindOneInt("seed", 0);


        return new OBDPTIntegrator(sampler, camera, minDepth, maxDepth, maxOptiDepth, pixelBounds, h, lightStrategy, conservative, seed);
    }

}  // namespace pbrt
