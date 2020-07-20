
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_BDPT_H
#define PBRT_INTEGRATORS_BDPT_H

// integrators/bdpt.h*
#include <unordered_map>
#include "camera.h"
#include "integrator.h"
#include "interaction.h"
#include "light.h"
#include "pbrt.h"
#include "reflection.h"
#include "sampling.h"
#include "scene.h"

// Moved shared bdpt functions to:
#include "bdpt_base.h"

namespace pbrt {

// BDPT Helper Definitions

struct Vertex;


// BDPT Declarations
class BDPTIntegrator : public Integrator {
  public:
    // BDPTIntegrator Public Methods
    BDPTIntegrator(std::shared_ptr<Sampler> sampler,
                   std::shared_ptr<const Camera> camera, int maxDepth,
                   bool visualizeStrategies, bool visualizeWeights,
                   const Bounds2i &pixelBounds,
                   const std::string &lightSampleStrategy = "power")
        : sampler(sampler),
          camera(camera),
          maxDepth(maxDepth),
          visualizeStrategies(visualizeStrategies),
          visualizeWeights(visualizeWeights),
          pixelBounds(pixelBounds),
          lightSampleStrategy(lightSampleStrategy) {}
    void Render(const Scene &scene);

  private:
    // BDPTIntegrator Private Data
    std::shared_ptr<Sampler> sampler;
    std::shared_ptr<const Camera> camera;
    const int maxDepth;
    const bool visualizeStrategies;
    const bool visualizeWeights;
    const Bounds2i pixelBounds;
    const std::string lightSampleStrategy;
};


Spectrum ConnectBDPT(
    const Scene &scene, Vertex *lightVertices, Vertex *cameraVertices, int s,
    int t, const Distribution1D &lightDistr,
    const std::unordered_map<const Light *, size_t> &lightToIndex,
    const Camera &camera, Sampler &sampler, Point2f *pRaster,
    Float *misWeight = nullptr);

BDPTIntegrator *CreateBDPTIntegrator(const ParamSet &params,
                                     std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera);



}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_BDPT_H
