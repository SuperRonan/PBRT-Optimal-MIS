
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

#ifndef PBRT_INTEGRATORS_OPTMIS_H
#define PBRT_INTEGRATORS_OPTMIS_H

// integrators/optmis.h*
#include "pbrt.h"
#include "integrator.h"
#include "scene.h"
#include "guiding.h"
#include "lightdistrib.h"
#include "Eigen/Dense"

namespace pbrt {

// MISWeights Declarations
enum class MISWeights { Balance, Power, Optimal };

// OptimalMode Declarations
// How the optimal weights are implemented
// AlphaSum = the result is the sum of estimated optimal alphas (corresponding to variance minimization)
// AlphaSumL2 = the result is the sum of estimated alphas corresponding to L2 minimization
// AlphaSumIndB = the result is the sum of estimated optimal alphas (corresponding to variance minimization), the A matrix and the b vector
//                use independent sets of samples 
enum class OptimalMode { AlphaSum, AlphaSumL2, AlphaSumIndB, FullWeights, Progressive, AlphaSumProg };

// SamplingTechnique Declarations
// Strategies for random direction selection
// B = brdf sampling, G = spherical light projection sampling, L = light surface sampling
enum class SamplingTechnique { B, G, L, BG, BL, GL, BGL };

// LightSelection Declarations
// Strategies for random light selection
// A = additional (defined by LightSelectAType), E = according to estimated contributions, U = uniform, ALL = no random selection, one sample from each light is taken
enum class LightSelection { A, E, U, AE, AU, EU, AEU, ALL };

// LightSelectAType Declarations
// Type of the A (additional) strategy for the random light selection
// Power = according to total power, Strongest = just the light with strongest estimated contribution, NotStrongest = according to estimated contributions ommiting the strongest one
enum class LightSelectAType { Power, Strongest, NotStrongest };

// AlphaEstimator Declarations
class AlphaEstimator
{
public:
	enum class Update { Matrix, RightSide, Both };
public:
	AlphaEstimator(int techniques, int updateStep, int skipTech) : 
		techniques(techniques), 
		updateStep(updateStep),
		skipTech(skipTech),
		decomp(techniques, techniques, Eigen::ComputeThinU | Eigen::ComputeThinV) {
		Reset();
	}

	void Reset() {
		A = Eigen::MatrixXf::Zero(techniques, techniques);
		b = Eigen::MatrixXf::Zero(techniques, 3);
		Aitc = Eigen::MatrixXf::Zero(techniques, 1);
		bitc = Eigen::MatrixXf::Zero(1, 3);

		alphas.clear();
		updates = 0;
	}

	void UpdateEstimates(const std::vector<Spectrum> &contribs, const std::vector<float> &pdfs, const std::vector<int> &sampleCounts, 
		MISWeights misWeights, OptimalMode optimalMode, Update update, bool debug);	

	Spectrum ComputeAlphaSum(bool debug) {
		if (alphas.empty()) {
			ComputeAlphas(debug);
		}		
		Spectrum sum(0.f);
		for (int i = 0; i < techniques + 1; i++) {
			sum += alphas[i];
		}
		return sum;
	}

	Spectrum GetAlpha(int i, bool debug) {
		if (alphas.empty()) {
			ComputeAlphas(debug);
		}
		return alphas[i];
	}

private:
	Eigen::MatrixXf A;
	Eigen::MatrixXf b;
	Eigen::MatrixXf Aitc;
	Eigen::MatrixXf bitc;
	Eigen::BDCSVD<Eigen::MatrixXf> decomp;
	std::vector<Spectrum> alphas;
	const int techniques;
	const int updateStep;
	const int skipTech;
	int updates;

	void ComputeAlphas(bool debug);
};

// OptimalMISIntegrator Declarations
class OptimalMISIntegrator : public Integrator {
public:
	struct Stats {
		std::vector<Spectrum> weights;
		std::vector<int> weightsUpdates;
		Spectrum sumOfContribsSquared;
		Spectrum sumOfContribs;
		int contribsCount;
        int missCount;
        int hitCount;

		Stats(const int techniques) {
			weights.clear();
			weightsUpdates.clear();
			for (int i = 0; i < techniques; i++) {
				weights.push_back(Spectrum(0.f));
				weightsUpdates.push_back(0);
			}
			
			sumOfContribs = Spectrum(0.f);
			sumOfContribsSquared = Spectrum(0.f);
			contribsCount = 0;
            missCount = 0;
            hitCount = 0;
		}
	};
public:
	// OptimalMISIntegrator Public Methods
	OptimalMISIntegrator(MISWeights misWeights, MISWeights misWeightsForAlpha, OptimalMode optimalMode, 
		SamplingTechnique samplingTechnique, LightSelection lightSelection, LightSelectAType lightSelectAType,
		GuidingDistribution::SamplingProjection lProjection,
		GuidingDistribution::SamplingProjection gProjection,
		int trainSamples,
		std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler, std::shared_ptr<Sampler> trainSampler,
		const Bounds2i &pixelBounds,
		int seedOffset, bool moreDataLayers, int updateStep, int skipTech, int targetTime)
		: misWeights(misWeights), misWeightsForAlpha(misWeightsForAlpha), optimalMode(optimalMode), 
		samplingTechnique(samplingTechnique), lightSelection(lightSelection), lightSelectAType(lightSelectAType),
		lProjection(lProjection), gProjection(gProjection), trainSamples(trainSamples), camera(camera),
		sampler(sampler), trainSampler(trainSampler), pixelBounds(pixelBounds), seedOffset(seedOffset), moreDataLayers(moreDataLayers), updateStep(updateStep), skipTech(skipTech),
		targetTime(targetTime) {}
	void Render(const Scene &scene);
	void RenderImpl(const Scene &scene, bool writeToImage, int64_t& renderTime);
	Spectrum Li(const RayDifferential &ray, const Scene &scene, Sampler &sampler, MemoryArena &arena, 
		std::vector<AlphaEstimator> &alphasVector, bool training, Stats *stats, bool debug) const;

private:
	// OptimalMISIntegrator Private Data
	const MISWeights misWeights;
	const MISWeights misWeightsForAlpha;
	const OptimalMode optimalMode;
	const SamplingTechnique samplingTechnique;
	LightSelection lightSelection;
	const LightSelectAType lightSelectAType;
	const GuidingDistribution::SamplingProjection lProjection; // projection used for the light surface sampling (L)
	const GuidingDistribution::SamplingProjection gProjection; // projection used for the spherical light projection sampling (G)
	const int trainSamples;
	std::shared_ptr<const Camera> camera;
	std::shared_ptr<Sampler> sampler;
	std::shared_ptr<Sampler> trainSampler;
	const Bounds2i pixelBounds;
	const int seedOffset;
	const bool moreDataLayers;
	const int updateStep;
	const int skipTech;
	const int targetTime;

	std::unique_ptr<LightDistribution> lightDistribUniform;
	std::unique_ptr<LightDistribution> lightDistribPower;
	std::unique_ptr<LightDistribution> lightDistribEstimates;

	// OptimalMISIntegrator Private Methods
	Spectrum EstimateLiOrAlphaCoefs(const SurfaceInteraction &it, const Light &light, const Scene &scene,
		const Point2f &uBSDF, const Point2f &uLight, const Point2f &uGuiding, AlphaEstimator& alphas, 
		float lightSelectPdf, float *pointSelectPd, bool training, Stats *stats, bool debug) const;

	Spectrum Weight(const std::vector<Spectrum> &contribs, const std::vector<float> &pdfs, 
		const std::vector<int> &sampleCounts, AlphaEstimator& alphas, Stats *stats, bool debug) const;
};

OptimalMISIntegrator *CreateOptimalMISIntegrator(
	const ParamSet &params, std::shared_ptr<Sampler> sampler,
	std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_OPTMIS_H
