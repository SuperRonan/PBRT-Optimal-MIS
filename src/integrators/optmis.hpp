
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


// integrators/optmis.cpp*
#include "integrators/optmis.h"
#include "interaction.h"
#include "paramset.h"
#include "progressreporter.h"
#include "camera.h"
#include "film.h"
#include "stats.h"
#include "samplers/random.h"
#include <chrono>
#include <fstream> 

namespace pbrt {

STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);

// AlphaEstimator Method Definitions
void AlphaEstimator::UpdateEstimates(const std::vector<Spectrum> &contribs, const std::vector<float> &pdfs, const std::vector<int> &sampleCounts, 
	MISWeights misWeights, OptimalMode optimalMode, Update update, bool debug) {
	if (debug) {
		LOG(ERROR) << "Alpha update";
		for (int i = 0; i < contribs.size(); ++i) {
			LOG(ERROR) << StringPrintf("contribs[%d]: ", i) << contribs[i].ToString() << "\n";
		}
		for (int i = 0; i < pdfs.size(); ++i) {
			LOG(ERROR) << StringPrintf("pdfs[%d]: %f\n", i, pdfs[i]);
		}

        LOG(ERROR) << "old A:\n" << A;
        LOG(ERROR) << "old b:\n" << b;
	}

	for (int i = 0; i < techniques; i++) {
		if (pdfs[i * techniques + i] == 0.f) return;
	}

	std::vector<float> weights(techniques);
	std::vector<float> weights2(techniques);
	for (int i = 0; i < techniques; i++) {
		float denom = 0.f;
		float denomSq = 0.f;
		for (int j = 0; j < techniques; j++) {
			denom += sampleCounts[j] * pdfs[i * techniques + j];
			denomSq += sampleCounts[j] * pdfs[i * techniques + j] * sampleCounts[j] * pdfs[i * techniques + j];
		}
		if (misWeights == MISWeights::Balance) {
			if (optimalMode == OptimalMode::AlphaSumL2) {
				weights[i] = sampleCounts[i] / denom;
			} else {
				weights[i] = sampleCounts[i] / (denom * denom);

				//for other options of misWeights/optimalMode, weights2 are undefined ..
				weights2[i] = 1 / denom;
			}
		} else {
			if (optimalMode == OptimalMode::AlphaSumL2) {
				weights[i] = sampleCounts[i] * sampleCounts[i] * pdfs[i * techniques + i] / denomSq;
			} else {
				weights[i] = sampleCounts[i] * sampleCounts[i] * pdfs[i * techniques + i] / (denomSq * denom);
			}
		}
	}

	for (int i = 0; i < techniques; i++) {
		if (update == Update::Both || update == Update::Matrix) {
			for (int j = 0; j < techniques; j++) {
				for (int k = 0; k < techniques; k++) {
					A(i, j) += weights[k] * pdfs[k * techniques + i] * pdfs[k * techniques + j];
				}
			}
		}
		if (update == Update::Both || update == Update::RightSide) {
			for (int c = 0; c < 3; ++c) {
				for (int k = 0; k < techniques; k++) {
					b(i, c) += weights[k] * float(contribs[k][c]) * pdfs[k * techniques + i];
				}
			}
		}
	}

	// intercept terms (Fan method)
	// defined only when (misWeights == MISWeights::Balance && optimalMode != OptimalMode::AlphaSumL2)
	if(skipTech != -1){
		for (int i = 0; i < techniques; i++) {
			if (update == Update::Both || update == Update::Matrix) {
				for (int k = 0; k < techniques; k++) {
					Aitc(i, 0) += weights2[k] * pdfs[k * techniques + i];
				}
			}
		}
		if (update == Update::Both || update == Update::RightSide) {
			for (int c = 0; c < 3; ++c) {
				for (int k = 0; k < techniques; k++) {
					bitc(0, c) += weights2[k] * float(contribs[k][c]);
				}
			}
		}
	}

	++updates;
	if (updates % updateStep == 0) {
		alphas.clear();
	}
	
    if (debug) {
        LOG(ERROR) << "new A:\n" << A;
        LOG(ERROR) << "new b:\n" << b;
    }
}

void AlphaEstimator::ComputeAlphas(bool debug) {
	if (debug) {
		LOG(ERROR) << "Alpha computation";
		LOG(ERROR) << "A:\n" << A;
		LOG(ERROR) << "b:\n" << b;
		LOG(ERROR) << "determinant A:" << A.determinant();
	}

	// accounts also for the intercept term
	alphas.resize(techniques + 1);

	if (updates == 0) {
		for (int i = 0; i < techniques + 1; i++) {
			alphas[i] = Spectrum(0.f);
		}
	}
	else {
		// our method
		if(skipTech == -1){
			decomp.compute(A);
			Eigen::MatrixXf alphasM = decomp.solve(b);
			
			for (int i = 0; i < techniques; i++) {
				for (int c = 0; c < 3; ++c) {
					alphas[i][c] = alphasM(i, c);
				}
			}

			if (debug) {
				LOG(ERROR) << "relative error:" << (A * alphasM - b).norm() / b.norm();
			}
		}
		// Fan et al. method
		else{
			// for debug
			// std::stringstream rr;

			// compose normal equation matrix along with the intercept terms
			Eigen::MatrixXf aa(techniques + 1, techniques + 1);
			Eigen::MatrixXf bb(techniques + 1, 3);
			// TODO: check the scale of the individual parts
			aa << 
				A * (techniques * techniques), Aitc * techniques, 
				Aitc.transpose() * techniques, updates * techniques;
			bb << b * techniques * techniques, bitc * techniques;

			// // for debug
			// {
			// 	rr << "matrix a:\n" << A << std::endl;
			// 	rr << "vector b:\n" << b << std::endl;

			// 	rr << "matrix aa:\n" << aa << std::endl;
			// 	rr << "vector bb:\n" << bb << std::endl;
			// 	//rr << "vector alphas:\n" << alphasM << std::endl;

				
			// 	rr << "matrix aa(middle):\n" << aa.middleRows(skipTech, techniques - skipTech) << std::endl;
			// 	rr << "matrix aa(bottom):\n" << aa.bottomRows(skipTech + 1) << std::endl;
			// }			
			
			// // skip one row and column
			aa.middleRows(skipTech, techniques - skipTech) = aa.bottomRows(techniques - skipTech).eval();
			aa.middleCols(skipTech, techniques - skipTech) = aa.rightCols(techniques - skipTech).eval();
			aa.conservativeResize(techniques, techniques);

			// skip one row
			bb.middleRows(skipTech, techniques - skipTech) = bb.bottomRows(techniques - skipTech).eval();
			bb.conservativeResize(techniques, 3);


			// // for debug
			// {
			// 	rr << "matrix aa(crop):\n" << aa << std::endl;
			// 	rr << "vector bb(crop):\n" << bb << std::endl;
			// 	//rr << "vector alphas:\n" << alphasM << std::endl;
			// }

			// solve
			decomp.compute(aa + Eigen::MatrixXf::Identity(techniques, techniques));
			// decomp.compute(aa + Eigen::MatrixXf::Identity(techniques+1, techniques+1));
			// decomp.compute(aa);
			// decomp.setThreshold(0.001f);
			Eigen::MatrixXf alphasM = decomp.solve(bb);
			
			

			for (int i = 0, ii = 0; ii < techniques + 1; i++, ii++) {
				if(ii == skipTech){
					alphas[ii] = Spectrum(0.f);
					i--;
				}
				else{
					for (int c = 0; c < 3; ++c) {
						alphas[ii][c] = alphasM(i, c);
					}
					
				}
			}
			
			// for (int i = 0, ii = 0; ii < techniques + 1; i++, ii++) {
			// 	for (int c = 0; c < 3; ++c) {
			// 		alphas[ii][c] = alphasM(i, c);
			// 	}
			// }

			// // for debug
			// LOG(ERROR) << rr.str();
		}		
	}

	if (debug) {
		for (int i = 0; i < techniques; i++) {
			LOG(ERROR) << StringPrintf("alpha[%d]: ", i) << alphas[i].ToString();
		}
	}
}

// OptimalMISIntegrator Method Definitions
Spectrum OptimalMISIntegrator::Li(const RayDifferential &ray, const Scene &scene, Sampler &sampler, 
	MemoryArena &arena, std::vector<AlphaEstimator> &alphasVector, bool training, Stats *stats, bool debug) const {
	ProfilePhase p(Prof::SamplerIntegratorLi);

	Spectrum L(0.f);
	// Find closest ray intersection or return background radiance
	SurfaceInteraction isect;
	if (!scene.Intersect(ray, &isect)) {
		for (const auto &light : scene.lights) L += light->Le(ray);
        stats->missCount++;
		return L;
	}

	// Compute scattering functions for surface interaction
	isect.ComputeScatteringFunctions(ray, arena);
	if (!isect.bsdf)
		return Li(isect.SpawnRay(ray.d), scene, sampler, arena, alphasVector, training, stats, debug);

    stats->hitCount++;

	// Compute emitted light if ray hit an area light source	
	const Vector3f wo = isect.wo;
	L += isect.Le(wo);

	// Sample direct illumination and update alpha estimates
	{
		ProfilePhase p(Prof::DirectLighting);

		int lightCount = int(scene.lights.size());
		if (lightCount == 1 || lightSelection == LightSelection::ALL) {
			float p_x;
			for (int l = 0; l < lightCount; ++l) {            
				L += EstimateLiOrAlphaCoefs(isect, *scene.lights[l], scene, sampler.Get2D(), sampler.Get2D(), sampler.Get2D(), alphasVector[l], 1.f, &p_x, training, stats, debug);
			}
		} else {
			const Distribution1D *distribU = lightDistribUniform->Lookup(isect.p);
			const Distribution1D *distribE = lightDistribEstimates->Lookup(isect.p);

			const Distribution1D *distribA = nullptr;
			std::unique_ptr<Distribution1D> distribAdata = nullptr;
			if (lightSelectAType == LightSelectAType::Power) {
				distribA = lightDistribPower->Lookup(isect.p);
			} else {
				std::vector<float> contribs(distribE->func);
				int maxContribIdx = 0;
				for (int i = 1; i < contribs.size(); ++i) {
					if (contribs[i] > contribs[maxContribIdx]) maxContribIdx = i;
				}
				if (lightSelectAType == LightSelectAType::NotStrongest) {
					contribs[maxContribIdx] = 0.f;
				} else {
					for (int i = 0; i < contribs.size(); ++i) {
						if (i != maxContribIdx) contribs[i] = 0.f;
					}
				}
				distribAdata = std::unique_ptr<Distribution1D>(new Distribution1D(&contribs[0], int(contribs.size())));
				distribA = distribAdata.get();
			}

            bool useA = (lightSelection == LightSelection::A || lightSelection == LightSelection::AE || lightSelection == LightSelection::AU || lightSelection == LightSelection::AEU);
            bool useE = (lightSelection == LightSelection::E || lightSelection == LightSelection::AE || lightSelection == LightSelection::EU || lightSelection == LightSelection::AEU);
            bool useU = (lightSelection == LightSelection::U || lightSelection == LightSelection::AU || lightSelection == LightSelection::EU || lightSelection == LightSelection::AEU);

			// 1) Select light uniformly

			Spectrum f_lU(0.f);
			float pU_lU = 0.f, pE_lU = 0.f, pA_lU = 0.f, p_xU = 0.f;
            if (useU) {
				const int lightIdx = distribU->SampleDiscrete(sampler.Get1D(), &pU_lU);
                if (useE) pE_lU = distribE->DiscretePDF(lightIdx);
                if (useA) pA_lU = distribA->DiscretePDF(lightIdx);
				f_lU = EstimateLiOrAlphaCoefs(isect, *scene.lights[lightIdx], scene, sampler.Get2D(), sampler.Get2D(), sampler.Get2D(), alphasVector[0], pU_lU, &p_xU, training, stats, debug);
			}			

			// 2) Select light according to its estimated contribution

			Spectrum f_lE(0.f);
			float pU_lE = 0.f, pE_lE = 0.f, pA_lE = 0.f, p_xE = 0.f;
			if (useE) {
				const int lightIdx = distribE->SampleDiscrete(sampler.Get1D(), &pE_lE);
                if (useU) pU_lE = distribU->DiscretePDF(lightIdx);
                if (useA) pA_lE = distribA->DiscretePDF(lightIdx);
				f_lE = EstimateLiOrAlphaCoefs(isect, *scene.lights[lightIdx], scene, sampler.Get2D(), sampler.Get2D(), sampler.Get2D(), alphasVector[0], pE_lE, &p_xE, training, stats, debug);
			}

			// 3) Select light according to its power / estimated contribution with the most important light excluded / just the most important light

			Spectrum f_lA(0.f);
			float pU_lA = 0.f, pE_lA = 0.f, pA_lA = 0.f, p_xA = 0.f;
			if (useA) {
				const int lightIdx = distribA->SampleDiscrete(sampler.Get1D(), &pA_lA);
                if (useU) pU_lA = distribU->DiscretePDF(lightIdx);
                if (useE) pE_lA = distribE->DiscretePDF(lightIdx);
				f_lA = EstimateLiOrAlphaCoefs(isect, *scene.lights[lightIdx], scene, sampler.Get2D(), sampler.Get2D(), sampler.Get2D(), alphasVector[0], pA_lA, &p_xA, training, stats, debug);
			}

			// Combine

			if (lightSelection == LightSelection::U) {
				L += f_lU;
			}
			else if (lightSelection == LightSelection::E) {
				L += f_lE;
			}
			else if (lightSelection == LightSelection::A) {
				L += f_lA;
			}
			else {
				std::vector<Spectrum> contribs;
				std::vector<float> pdfs;
				std::vector<int> sampleCounts;

				switch (lightSelection) {
				case LightSelection::AE:
					contribs = { f_lA, f_lE };
					pdfs = { pA_lA * p_xA, pE_lA * p_xA, pA_lE * p_xE, pE_lE * p_xE };
					sampleCounts = { 1, 1 };
					break;
				case LightSelection::AU:
					contribs = { f_lA, f_lU };
					pdfs = { pA_lA * p_xA, pU_lA * p_xA, pA_lU * p_xU, pU_lU * p_xU };
					sampleCounts = { 1, 1 };
					break;
				case LightSelection::EU:
					contribs = { f_lE, f_lU };
					pdfs = { pE_lE * p_xE, pU_lE * p_xE, pE_lU * p_xU, pU_lU * p_xU };
					sampleCounts = { 1, 1 };
					break;
				case LightSelection::AEU:
					contribs = { f_lA, f_lE, f_lU };
					pdfs = { pA_lA * p_xA, pE_lA * p_xA, pU_lA * p_xA, pA_lE * p_xE, pE_lE * p_xE, pU_lE * p_xE, pA_lU * p_xU, pE_lU * p_xU, pU_lU * p_xU };
					sampleCounts = { 1, 1, 1 };
					break;
				}

                if (misWeights == MISWeights::Optimal && optimalMode == OptimalMode::Progressive) {
                    L += Weight(contribs, pdfs, sampleCounts, alphasVector[0], stats, debug);
                    alphasVector[0].UpdateEstimates(contribs, pdfs, sampleCounts, misWeightsForAlpha, optimalMode, AlphaEstimator::Update::Both, debug);
                } else {
                    if (misWeights != MISWeights::Optimal || (optimalMode == OptimalMode::FullWeights && !training)) {
                        L += Weight(contribs, pdfs, sampleCounts, alphasVector[0], stats, debug);
                    } else {
                        if (optimalMode != OptimalMode::AlphaSumIndB) {
                            alphasVector[0].UpdateEstimates(contribs, pdfs, sampleCounts, misWeightsForAlpha, optimalMode, AlphaEstimator::Update::Both, debug);
                        } else {
                            // If the matrix and the right side should use independent samples, update first the matrix with the current samples
                            alphasVector[0].UpdateEstimates(contribs, pdfs, sampleCounts, misWeightsForAlpha, optimalMode, AlphaEstimator::Update::Matrix, debug);

                            // Generate second set of samples
                            f_lU = Spectrum(0.f);
                            pU_lU = 0.f, pE_lU = 0.f, pA_lU = 0.f, p_xU = 0.f;
                            if (lightSelection == LightSelection::U || lightSelection == LightSelection::AU || lightSelection == LightSelection::EU || lightSelection == LightSelection::AEU) {
                                const int lightIdx = distribU->SampleDiscrete(sampler.Get1D(), &pU_lU);
                                pE_lU = distribE->DiscretePDF(lightIdx);
                                pA_lU = distribA->DiscretePDF(lightIdx);
                                f_lU = EstimateLiOrAlphaCoefs(isect, *scene.lights[lightIdx], scene, sampler.Get2D(), sampler.Get2D(), sampler.Get2D(), alphasVector[0], pU_lU, &p_xU, training, stats, debug);
                            }
                            f_lE = Spectrum(0.f);
                            pU_lE = 0.f, pE_lE = 0.f, pA_lE = 0.f, p_xE = 0.f;
                            if (lightSelection == LightSelection::E || lightSelection == LightSelection::AE || lightSelection == LightSelection::EU || lightSelection == LightSelection::AEU) {
                                const int lightIdx = distribE->SampleDiscrete(sampler.Get1D(), &pE_lE);
                                pU_lE = distribU->DiscretePDF(lightIdx);
                                pA_lE = distribA->DiscretePDF(lightIdx);
                                f_lE = EstimateLiOrAlphaCoefs(isect, *scene.lights[lightIdx], scene, sampler.Get2D(), sampler.Get2D(), sampler.Get2D(), alphasVector[0], pE_lE, &p_xE, training, stats, debug);
                            }
                            f_lA = Spectrum(0.f);
                            pU_lA = 0.f, pE_lA = 0.f, pA_lA = 0.f, p_xA = 0.f;
                            if (lightSelection == LightSelection::A || lightSelection == LightSelection::AE || lightSelection == LightSelection::AU || lightSelection == LightSelection::AEU) {
                                const int lightIdx = distribA->SampleDiscrete(sampler.Get1D(), &pA_lA);
                                pU_lA = distribU->DiscretePDF(lightIdx);
                                pE_lA = distribE->DiscretePDF(lightIdx);
                                f_lA = EstimateLiOrAlphaCoefs(isect, *scene.lights[lightIdx], scene, sampler.Get2D(), sampler.Get2D(), sampler.Get2D(), alphasVector[0], pA_lA, &p_xA, training, stats, debug);
                            }
                            switch (lightSelection) {
                            case LightSelection::AE:
                                contribs = { f_lA, f_lE };
                                pdfs = { pA_lA * p_xA, pE_lA * p_xA, pA_lE * p_xE, pE_lE * p_xE };
                                sampleCounts = { 1, 1 };
                                break;
                            case LightSelection::AU:
                                contribs = { f_lA, f_lU };
                                pdfs = { pA_lA * p_xA, pU_lA * p_xA, pA_lU * p_xU, pU_lU * p_xU };
                                sampleCounts = { 1, 1 };
                                break;
                            case LightSelection::EU:
                                contribs = { f_lE, f_lU };
                                pdfs = { pE_lE * p_xE, pU_lE * p_xE, pE_lU * p_xU, pU_lU * p_xU };
                                sampleCounts = { 1, 1 };
                                break;
                            case LightSelection::AEU:
                                contribs = { f_lA, f_lE, f_lU };
                                pdfs = { pA_lA * p_xA, pE_lA * p_xA, pU_lA * p_xA, pA_lE * p_xE, pE_lE * p_xE, pU_lE * p_xE, pA_lU * p_xU, pE_lU * p_xU, pU_lU * p_xU };
                                sampleCounts = { 1, 1, 1 };
                                break;
                            }

                            // Update the right side with the second set
                            alphasVector[0].UpdateEstimates(contribs, pdfs, sampleCounts, misWeightsForAlpha, optimalMode, AlphaEstimator::Update::RightSide, debug);
                        }
                    }
                }
            }
        }
	}

	return L;
}

Spectrum OptimalMISIntegrator::EstimateLiOrAlphaCoefs(const SurfaceInteraction &it, const Light &light, const Scene &scene, 
	const Point2f &uBSDF, const Point2f &uLight, const Point2f &uGuiding, AlphaEstimator& alphas, 
	float lightSelectPdf, float *pointSelectPdf, bool training, Stats *stats, bool debug) const {
	// todo: delta lights, samples with zero pdfs

	BxDFType bsdfFlags = BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
	GuidingDistribution guiding(it, light);
	bool normalLightSampling = (lProjection == GuidingDistribution::SamplingProjection::None) || !guiding.CanSample(lProjection);
	Vector3f wi;
	Spectrum Li;
	VisibilityTester visibility;
    
    bool useB = (samplingTechnique == SamplingTechnique::B || samplingTechnique == SamplingTechnique::BL || samplingTechnique == SamplingTechnique::BG || samplingTechnique == SamplingTechnique::BGL);
    bool useG = (samplingTechnique == SamplingTechnique::G || samplingTechnique == SamplingTechnique::BG || samplingTechnique == SamplingTechnique::GL || samplingTechnique == SamplingTechnique::BGL);
    bool useL = (samplingTechnique == SamplingTechnique::L || samplingTechnique == SamplingTechnique::BL || samplingTechnique == SamplingTechnique::GL || samplingTechnique == SamplingTechnique::BGL);

	// 1) Light sampling

	Spectrum f_xL(0.f);
	Float pL_xL = 0, pB_xL = 0, pG_xL = 0;

    if (useL) {
        if (normalLightSampling) {
            // Sample original light -> Li, pL_xL
            Li = light.Sample_Li(it, uLight, &wi, &pL_xL, &visibility);

            // If the light sample is valid
            if (pL_xL > 0.f) {
                // Compute pdf of BSDF sampling -> pB_xL
                if (useB) pB_xL = it.bsdf->Pdf(it.wo, wi, bsdfFlags);

                // Compute pdf of guiding -> pG_xL
                if (useG) pG_xL = guiding.Pdf(wi, gProjection);

                if (!Li.IsBlack()) {
                    // Compute BSDF -> f*dot
                    Spectrum f = it.bsdf->f(it.wo, wi, bsdfFlags) * AbsDot(wi, it.shading.n);

                    // Compute sample contribution -> f_xL
                    if (!f.IsBlack() && visibility.Unoccluded(scene)) {
                        f_xL = f * Li;
                    }
                }
            }
        } else {
            // Sampling projected light -> pL_xL
            wi = guiding.Sample_wi(uLight, lProjection, &pL_xL);

            // If the light sample is valid
            if (pL_xL > 0.f) {
                // Compute pdf of BSDF sampling -> pB_xG
                if (useB) pB_xL = it.bsdf->Pdf(it.wo, wi, bsdfFlags);

                // Compute pdf of guiding -> pG_xL
                if (useG) pG_xL = guiding.Pdf(wi, gProjection);

                // Compute BSDF -> f*dot
                Spectrum f = it.bsdf->f(it.wo, wi, bsdfFlags) * AbsDot(wi, it.shading.n);

                if (!f.IsBlack()) {
                    // Compute incoming radiance -> Li
                    SurfaceInteraction lightIsect;
                    Ray ray = it.SpawnRay(wi);
                    Li = Spectrum(0.f);
                    if (scene.Intersect(ray, &lightIsect)) {
                        if (lightIsect.primitive->GetAreaLight() == &light) {
                            Li = lightIsect.Le(-wi);
                        }
                    } else {
                        Li = light.Le(ray);
                    }

                    // Compute sample contribution -> f_xG
                    f_xL = f * Li;
                }
            }
        }
    }

	// 2) BSDF sampling

	Spectrum f_xB(0.f);
	Float pL_xB = 0, pB_xB = 0, pG_xB = 0;

    if (useB) {
        // Sample BSDF -> f*dot, pB_xB
        BxDFType sampledType;
        Spectrum f = it.bsdf->Sample_f(it.wo, &wi, uBSDF, &pB_xB, bsdfFlags, &sampledType);
        f *= AbsDot(wi, it.shading.n);

        // If the BSDF sample is valid
        if (pB_xB > 0.f) {
            // Compute pdf of light sampling -> pL_xB
            if (useL) {
                if (normalLightSampling) {
                    pL_xB = light.Pdf_Li(it, wi);
                } else {
                    pL_xB = guiding.Pdf(wi, lProjection);
                }
            }

            // Compute pdf of guiding -> pG_xL
            if (useG) pG_xB = guiding.Pdf(wi, gProjection);

            if (!f.IsBlack()) {
                // Compute incoming radiance -> Li
                SurfaceInteraction lightIsect;
                Ray ray = it.SpawnRay(wi);
                Li = Spectrum(0.f);
                if (scene.Intersect(ray, &lightIsect)) {
                    if (lightIsect.primitive->GetAreaLight() == &light) {
                        Li = lightIsect.Le(-wi);
                    }
                } else {
                    Li = light.Le(ray);
                }

                // Compute sample contribution -> f_xB
                f_xB = f * Li;
            }
        }
    }

	// 3) Guiding distribution sampling

	Spectrum f_xG(0.f);
	Float pL_xG = 0, pB_xG = 0, pG_xG = 0;

    if (useG) {
        // Sample guiding distribution -> wi, pG_xG
        wi = guiding.Sample_wi(uGuiding, gProjection, &pG_xG);

        // If the guiding sample is valid
        if (pG_xG > 0.f) {
            // Compute pdf of light sampling -> pL_xG
            if (useL) {
                if (normalLightSampling) {
                    pL_xG = light.Pdf_Li(it, wi);
                } else {
                    pL_xG = guiding.Pdf(wi, lProjection);
                }
            }

            // Compute pdf of BSDF sampling -> pB_xG
            if (useB) pB_xG = it.bsdf->Pdf(it.wo, wi, bsdfFlags);

            // Compute BSDF -> f*dot
            Spectrum f = it.bsdf->f(it.wo, wi, bsdfFlags) * AbsDot(wi, it.shading.n);

            if (!f.IsBlack()) {
                // Compute incoming radiance -> Li
                SurfaceInteraction lightIsect;
                Ray ray = it.SpawnRay(wi);
                Li = Spectrum(0.f);
                if (scene.Intersect(ray, &lightIsect)) {
                    if (lightIsect.primitive->GetAreaLight() == &light) {
                        Li = lightIsect.Le(-wi);
                    }
                } else {
                    Li = light.Le(ray);
                }

                // Compute sample contribution -> f_xG
                f_xG = f * Li;
            }
        }
    }

	// Combine

	if (lightSelection != LightSelection::A && lightSelection != LightSelection::E && lightSelection != LightSelection::U && lightSelection != LightSelection::ALL) {
		switch (samplingTechnique) {
		case SamplingTechnique::B:
			*pointSelectPdf = pB_xB;
			return f_xB;
		case SamplingTechnique::G:
			*pointSelectPdf = pG_xG;
			return f_xG;
		case SamplingTechnique::L:
			*pointSelectPdf = pL_xL;
			return f_xL;
		}
	} else {
		std::vector<Spectrum> contribs;
		std::vector<float> pdfs;
		std::vector<int> sampleCounts;

		switch (samplingTechnique) {
		case SamplingTechnique::B:
			contribs = { f_xB };
			pdfs = { pB_xB };
			sampleCounts = { 1 };
			break;
		case SamplingTechnique::G:
			contribs = { f_xG };
			pdfs = { pG_xG };
			sampleCounts = { 1 };
			break;
		case SamplingTechnique::L:
			contribs = { f_xL };
			pdfs = { pL_xL };
			sampleCounts = { 1 };
			break;
		case SamplingTechnique::BG:
			contribs = { f_xG, f_xB };
			pdfs = { pG_xG, pB_xG, pG_xB, pB_xB };
			sampleCounts = { 1, 1 };
			break;
		case SamplingTechnique::BL:
			contribs = { f_xL, f_xB };
			pdfs = { pL_xL, pB_xL, pL_xB, pB_xB };
			sampleCounts = { 1, 1 };
			break;
		case SamplingTechnique::GL:
			contribs = { f_xL, f_xG };
			pdfs = { pL_xL, pG_xL, pL_xG, pG_xG };
			sampleCounts = { 1, 1 };
			break;
		case SamplingTechnique::BGL:
			contribs = { f_xL, f_xB, f_xG };
			pdfs = { pL_xL, pB_xL, pG_xL, pL_xB, pB_xB, pG_xB, pL_xG, pB_xG, pG_xG };
			sampleCounts = { 1, 1, 1 };
			break;
		}

		for (int i = 0; i < pdfs.size(); ++i) {
			pdfs[i] *= lightSelectPdf;
		}

        if (misWeights == MISWeights::Optimal && optimalMode == OptimalMode::Progressive) {
            Spectrum L = Weight(contribs, pdfs, sampleCounts, alphas, stats, debug);
            alphas.UpdateEstimates(contribs, pdfs, sampleCounts, misWeightsForAlpha, optimalMode, AlphaEstimator::Update::Both, debug);
            return L;
        } else {
            if (misWeights != MISWeights::Optimal || (optimalMode == OptimalMode::FullWeights && !training)) {
                return Weight(contribs, pdfs, sampleCounts, alphas, stats, debug);
            } else {
                alphas.UpdateEstimates(contribs, pdfs, sampleCounts, misWeightsForAlpha, optimalMode, AlphaEstimator::Update::Both, debug);
                return Spectrum(0.f);
            }
        }
	}
}

OptimalMISIntegrator *CreateOptimalMISIntegrator(
	const ParamSet &params, std::shared_ptr<Sampler> sampler,
	std::shared_ptr<const Camera> camera) {
	MISWeights misWeights;
	std::string st = params.FindOneString("misweights", "optimal");
	if (st == "balance")
		misWeights = MISWeights::Balance;
	else if (st == "power")
		misWeights = MISWeights::Power;
	else if (st == "optimal")
		misWeights = MISWeights::Optimal;
	else {
		Warning("MIS weights \"%s\" unknown, using \"optimal\".", st.c_str());
		misWeights = MISWeights::Optimal;
	}
	MISWeights misWeightsForAlpha;
	st = params.FindOneString("misweightsa", "balance");
	if (st == "balance")
		misWeightsForAlpha = MISWeights::Balance;
	else if (st == "power")
		misWeightsForAlpha = MISWeights::Power;
	else {
		Warning("MIS weights for alpha \"%s\" unknown, using \"balance\".", st.c_str());
		misWeightsForAlpha = MISWeights::Optimal;
	}
	OptimalMode optimalMode;
	st = params.FindOneString("optimalmode", "alphasum");
	if (st == "alphasum")
		optimalMode = OptimalMode::AlphaSum;
	else if (st == "alphasuml2")
		optimalMode = OptimalMode::AlphaSumL2;
	else if (st == "alphasumindb")
		optimalMode = OptimalMode::AlphaSumIndB; 
	else if (st == "fullweights")
		optimalMode = OptimalMode::FullWeights;
    else if (st == "progressive")
        optimalMode = OptimalMode::Progressive;
    else if (st == "alphasumprog")
        optimalMode = OptimalMode::AlphaSumProg;
	else {
		Warning("Optimal mode \"%s\" unknown, using \"alphasum\".", st.c_str());
		optimalMode = OptimalMode::AlphaSum;
	}
	SamplingTechnique samplingTechnique;
	st = params.FindOneString("technique", "BGL");
	if (st == "B")
		samplingTechnique = SamplingTechnique::B;
	else if (st == "G")
		samplingTechnique = SamplingTechnique::G;
	else if (st == "L")
		samplingTechnique = SamplingTechnique::L;
	else if (st == "BG")
		samplingTechnique = SamplingTechnique::BG;
	else if (st == "BL")
		samplingTechnique = SamplingTechnique::BL;
	else if (st == "GL")
		samplingTechnique = SamplingTechnique::GL;
	else if (st == "BGL")
		samplingTechnique = SamplingTechnique::BGL;
	else {
		Warning("Sampling technique \"%s\" unknown, using \"BGL\".", st.c_str());
		samplingTechnique = SamplingTechnique::BGL;
	}
	LightSelection lightSelection;
	st = params.FindOneString("lselect", "U");
	if (st == "A")
		lightSelection = LightSelection::A;
	else if (st == "E")
		lightSelection = LightSelection::E;
	else if (st == "U")
		lightSelection = LightSelection::U;
	else if (st == "AE")
		lightSelection = LightSelection::AE;
	else if (st == "AU")
		lightSelection = LightSelection::AU;
	else if (st == "EU")
		lightSelection = LightSelection::EU;
	else if (st == "AEU")
		lightSelection = LightSelection::AEU;
	else if (st == "ALL")
		lightSelection = LightSelection::ALL;
	else {
		Warning("Light selection \"%s\" unknown, using \"U\".", st.c_str());
		lightSelection = LightSelection::U;
	}
	LightSelectAType lightSelectAType;
	st = params.FindOneString("lsat", "NS");
	if (st == "P")
		lightSelectAType = LightSelectAType::Power;
	else if (st == "S")
		lightSelectAType = LightSelectAType::Strongest;
	else if (st == "NS")
		lightSelectAType = LightSelectAType::NotStrongest;
	else {
		Warning("Light selection strategy A \"%s\" unknown, using \"NS\".", st.c_str());
		lightSelectAType = LightSelectAType::NotStrongest;
	}
	GuidingDistribution::SamplingProjection lProjection;
	st = params.FindOneString("lproject", "none");
	if (st == "none")
		lProjection = GuidingDistribution::SamplingProjection::None;
	else if (st == "parallelplane")
		lProjection = GuidingDistribution::SamplingProjection::ParallelPlane;
	else {
		Warning("Light sampling projection \"%s\" unknown, using \"none\".", st.c_str());
		lProjection = GuidingDistribution::SamplingProjection::None;
	}
	GuidingDistribution::SamplingProjection gProjection;
	st = params.FindOneString("gproject", "sphereprecise");
	if (st == "sphereprecise")
		gProjection = GuidingDistribution::SamplingProjection::SpherePrecise;
	else if (st == "spheresimple")
		gProjection = GuidingDistribution::SamplingProjection::SphereSimple;
	else {
		Warning("Guided sampling projection \"%s\" unknown, using \"sphereprecise\".", st.c_str());
		gProjection = GuidingDistribution::SamplingProjection::SpherePrecise;
	}
	int trainSamples = params.FindOneInt("trainsamples", 0);
	std::shared_ptr<Sampler> trainSampler = trainSamples == 0 ? nullptr : std::shared_ptr<Sampler>(new RandomSampler(trainSamples));
	int np;
	const int *pb = params.FindInt("pixelbounds", &np);
	Bounds2i pixelBounds = camera->film->GetSampleBounds();
	if (pb) {
		if (np != 4)
			Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
			np);
		else {
			pixelBounds = Intersect(pixelBounds,
				Bounds2i{ { pb[0], pb[2] }, { pb[1], pb[3] } });
			if (pixelBounds.Area() == 0)
				Error("Degenerate \"pixelbounds\" specified.");
		}
	}
	int seedOffset = params.FindOneInt("seedoffset", 0);
	bool moreDataLayers = params.FindOneBool("moredatalayers", false);
	int updateStep = params.FindOneInt("updatestep", 1);
	int skipTech = params.FindOneInt("skipTech", -1);
	int targetTime = params.FindOneInt("targettime", -1);
	return new OptimalMISIntegrator(misWeights, misWeightsForAlpha, optimalMode, samplingTechnique, lightSelection, lightSelectAType,
		lProjection, gProjection, trainSamples, camera, sampler, trainSampler, pixelBounds, seedOffset, moreDataLayers, updateStep, skipTech, targetTime);
}

void OptimalMISIntegrator::Render(const Scene &scene) {
	if (targetTime < 1) {
		int64_t renderTime;
		RenderImpl(scene, true, renderTime);
	} else {        
		int64_t renderTime;
		sampler->samplesPerPixel = 10;
		RenderImpl(scene, false, renderTime);

		sampler->samplesPerPixel = (int64_t)std::round(10 * (targetTime / (0.001f * ((Float)renderTime))));
		RenderImpl(scene, true, renderTime);
	}
}

void OptimalMISIntegrator::RenderImpl(const Scene &scene, bool writeToImage, int64_t& renderTime) {
    /*std::ofstream statsFileTimes;
    statsFileTimes.open("statsFileTimes.txt", std::ios::app);

    std::ofstream statsFileSamples;
    statsFileSamples.open("statsFileSamples.txt", std::ios::app);*/

    
    // Render image tiles in parallel

	// Compute number of tiles, _nTiles_, to use for parallel rendering
	Bounds2i sampleBounds = camera->film->GetSampleBounds();
	Vector2i sampleExtent = sampleBounds.Diagonal();
	const int tileSize = 16;
	Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
		(sampleExtent.y + tileSize - 1) / tileSize);

	if (scene.lights.size() >= 1) {
		if (scene.lights.size() == 1) {
			lightSelection = LightSelection::U;
		}
		else {
			if ((lightSelection != LightSelection::A && lightSelection != LightSelection::E && lightSelection != LightSelection::U && lightSelection != LightSelection::ALL) &&
				(samplingTechnique != SamplingTechnique::B && samplingTechnique != SamplingTechnique::G && samplingTechnique != SamplingTechnique::L)) {
				LOG(ERROR) << "Cannot combine more light selection techniques with more light sampling techniques. Setting light selection to U.";
				lightSelection = LightSelection::U;
			}
		}

		// Prepare light selection distributions
		lightDistribUniform = CreateLightSampleDistribution("uniform", scene);
		lightDistribPower = CreateLightSampleDistribution("power", scene);
		lightDistribEstimates = CreateLightSampleDistribution("spatial", scene);

		// Number of techniques
		int techniques;
		if (lightSelection == LightSelection::A || lightSelection == LightSelection::E || lightSelection == LightSelection::U || lightSelection == LightSelection::ALL) {
			if (samplingTechnique == SamplingTechnique::B || samplingTechnique == SamplingTechnique::G || samplingTechnique == SamplingTechnique::L) {
				techniques = 1;
			}
			else if (samplingTechnique != SamplingTechnique::BGL) {
				techniques = 2;
			}
			else {
				techniques = 3;
			}
		}
		else if (lightSelection != LightSelection::AEU) {
			techniques = 2;
		}
		else {
			techniques = 3;
		}

		// Extra exr data channels
		if (moreDataLayers && writeToImage) {
			// Weights
			for (int i = 0; i < techniques; ++i) {
				camera->film->dataNames.push_back("W" + std::to_string(i) + ".R");
				camera->film->dataNames.push_back("W" + std::to_string(i) + ".G");
				camera->film->dataNames.push_back("W" + std::to_string(i) + ".B");
			}
			// Alphas
			for (int i = 0; i < techniques; ++i) {
				camera->film->dataNames.push_back("A" + std::to_string(i) + ".R");
				camera->film->dataNames.push_back("A" + std::to_string(i) + ".G");
				camera->film->dataNames.push_back("A" + std::to_string(i) + ".B");
			}
			// Variance
			camera->film->dataNames.push_back("V.R");
			camera->film->dataNames.push_back("V.G");
			camera->film->dataNames.push_back("V.B");
		}

		const std::chrono::system_clock::time_point renderingStart = std::chrono::system_clock::now();

		ProgressReporter reporter(nTiles.x * nTiles.y, writeToImage ? "Rendering" : "Dry run");
		{
			ParallelFor2D([&](Point2i tile) {
				// Render section of image corresponding to _tile_

				// Allocate _MemoryArena_ for tile
				MemoryArena arena;

				// Get sampler instance for tile
				int seed = tile.y * nTiles.x + tile.x + seedOffset;
				std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);
				std::unique_ptr<Sampler> trainTileSampler = trainSamples == 0 ? nullptr : trainSampler->Clone(seed + 721);

				// Compute sample bounds for tile
				int x0 = sampleBounds.pMin.x + tile.x * tileSize;
				int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
				int y0 = sampleBounds.pMin.y + tile.y * tileSize;
				int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
				Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
				LOG(INFO) << "Starting image tile " << tileBounds;

				// Get _FilmTile_ for tile
				std::unique_ptr<FilmTile> filmTile =
					camera->film->GetFilmTile(tileBounds);

                std::vector<AlphaEstimator> alphasVector;
                if (lightSelection == LightSelection::ALL) {
                    for (int l = 0; l < scene.lights.size(); ++l) alphasVector.emplace_back(techniques, updateStep, skipTech);
                } else {
                    alphasVector.emplace_back(techniques, updateStep, skipTech);
                }
				
                // Loop over pixels in tile to render them
				for (Point2i pixel : tileBounds) {
					{
						ProfilePhase pp(Prof::StartPixel);
                        if (optimalMode == OptimalMode::Progressive || optimalMode == OptimalMode::AlphaSumProg) {
                            tileSampler = sampler->Clone(seed + pixel.x + pixel.y * camera->film->fullResolution.x);
                        }
						tileSampler->StartPixel(pixel);
						if (trainTileSampler) trainTileSampler->StartPixel(pixel);
					}

					bool debug = false;
					/*if ((pixel.x == 326 || pixel.x == 327) && pixel.y == 521) {
						debug = true;
						LOG(ERROR) << StringPrintf("pixel [%d, %d]\n", pixel.x, pixel.y);
					}*/

					// Do this check after the StartPixel() call; this keeps
					// the usage of RNG values from (most) Samplers that use
					// RNGs consistent, which improves reproducability /
					// debugging.
					if (!InsideExclusive(pixel, pixelBounds))
						continue;

                    for (int l = 0; l < alphasVector.size(); l++) alphasVector[l].Reset();

					Stats stats(techniques);

					//
					// 1. Perform alpha estimation on the training samples
					//

					if (trainSamples > 0) {
						do {
							// Initialize _CameraSample_ for current sample
							CameraSample cameraSample =
								trainTileSampler->GetCameraSample(pixel);
							

							if (debug) {
								LOG(ERROR) << StringPrintf("training sample %d\n", (int)trainTileSampler->CurrentSampleNumber());
							}

							// Dirty hack to avoid one alpha sum being added to two adjacent pixels if their exact boundary is sampled
							if (cameraSample.pFilm[0] == pixel.x) {
								cameraSample.pFilm[0] += 0.0001;
							}
							else if (cameraSample.pFilm[0] == pixel.x + 1) {
								cameraSample.pFilm[0] -= 0.0001;
							}
							if (cameraSample.pFilm[1] == pixel.y) {
								cameraSample.pFilm[1] += 0.0001;
							}
							else if (cameraSample.pFilm[1] == pixel.y + 1) {
								cameraSample.pFilm[1] -= 0.0001;
							}

							// Generate camera ray for current sample
							RayDifferential ray;
							Float rayWeight =
								camera->GenerateRayDifferential(cameraSample, &ray);

							// Hack to have always the same differentials scaling (to ease comparisons with reference etc ..)
							/*ray.ScaleDifferentials(
								1 / std::sqrt((Float)trainTileSampler->samplesPerPixel));*/
							ray.ScaleDifferentials(1 / std::sqrt(20));
							++nCameraRays;

							// Update pixel alpha estimates
							if (rayWeight > 0) Li(ray, scene, *trainTileSampler, arena, alphasVector, true, &stats, debug);

							// Free _MemoryArena_ memory from computing image sample
							// value
							arena.Reset();
						} while (trainTileSampler->StartNextSample());
					}

					//
					// 2. Perform DI/alpha estimation on the target samples
					//
					
					// Initialize _CameraSample_ for all samples

					// If Fan et al., do it just once here, so that several iterations count as one theirs. 
					// This hack is due to our inability to specify a different number of samples per technique than
					// one.
					CameraSample cameraSample =
						tileSampler->GetCameraSample(pixel);

					do {
						// If not Fan et. al .. Initialize _CameraSample_ always
						if(skipTech == -1)
							cameraSample =
								tileSampler->GetCameraSample(pixel);

						if (debug) {
							LOG(ERROR) << StringPrintf("sample %d\n", (int)tileSampler->CurrentSampleNumber());
						}

						// Dirty hack to avoid one alpha sum being added to two adjacent pixels if their exact boundary is sampled
						if (cameraSample.pFilm[0] == pixel.x) {
							cameraSample.pFilm[0] += 0.0001;
						}
						else if (cameraSample.pFilm[0] == pixel.x + 1) {
							cameraSample.pFilm[0] -= 0.0001;
						}
						if (cameraSample.pFilm[1] == pixel.y) {
							cameraSample.pFilm[1] += 0.0001;
						}
						else if (cameraSample.pFilm[1] == pixel.y + 1) {
							cameraSample.pFilm[1] -= 0.0001;
						}

						// Generate camera ray for current sample
						RayDifferential ray;
						Float rayWeight =
							camera->GenerateRayDifferential(cameraSample, &ray);

						// Hack to have always the same differentials scaling (to ease comparisons with reference etc ..)
						/*ray.ScaleDifferentials(
							1 / std::sqrt((Float)tileSampler->samplesPerPixel));*/
						ray.ScaleDifferentials(1 / std::sqrt(20));
							
						++nCameraRays;

						// Evaluate directly visible radiance along camera ray and update pixel alpha estimates
						Spectrum L(0.f);
						if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena, alphasVector, false, &stats, debug);

						// Issue warning if unexpected radiance value returned
						if (L.HasNaNs()) {
							LOG(ERROR) << StringPrintf(
								"Not-a-number radiance value returned "
								"for pixel (%d, %d), sample %d. Setting to black.",
								pixel.x, pixel.y,
								(int)tileSampler->CurrentSampleNumber());
							L = Spectrum(0.f);
						}
						else if (L.y() < -1e-5) {
							//LOG(ERROR) << StringPrintf(
							//	"Negative luminance value, %f, returned "
							//	"for pixel (%d, %d), sample %d. Setting to black.",
							//	L.y(), pixel.x, pixel.y,
							//	(int)tileSampler->CurrentSampleNumber());
							//L = Spectrum(0.f);
						}
						else if (std::isinf(L.y())) {
							LOG(ERROR) << StringPrintf(
								"Infinite luminance value returned "
								"for pixel (%d, %d), sample %d. Setting to black.",
								pixel.x, pixel.y,
								(int)tileSampler->CurrentSampleNumber());
							L = Spectrum(0.f);
						}
						VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
							ray << " -> L = " << L;

						// Add camera ray's contribution to image
                        filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

						// Free _MemoryArena_ memory from computing image sample
						// value
						arena.Reset();
					} while (tileSampler->StartNextSample());

					//
					// 3. Compute the sum of alphas and add it to the resulting image
					//

                    if (misWeights == MISWeights::Optimal && optimalMode != OptimalMode::FullWeights && optimalMode != OptimalMode::Progressive)
					{
						// Evaluate direct illumination for the current pixel
                        Spectrum alphaSum(0.f);
                        for (int l = 0; l < alphasVector.size(); l++)  alphaSum += alphasVector[l].ComputeAlphaSum(debug);

						// Issue warning if unexpected radiance value returned
						if (alphaSum.HasNaNs()) {
							LOG(ERROR) << StringPrintf(
								"Not-a-number alpha sum returned "
								"for pixel (%d, %d), sample %d. Setting to black.",
								pixel.x, pixel.y,
								(int)tileSampler->CurrentSampleNumber());
							alphaSum = Spectrum(0.f);
						}
						else if (alphaSum.y() < -1e-5) {
							LOG(ERROR) << StringPrintf(
								"Negative alpha sum, %f, returned "
								"for pixel (%d, %d), sample %d. Setting to black.",
								alphaSum.y(), pixel.x, pixel.y,
								(int)tileSampler->CurrentSampleNumber());
							alphaSum = Spectrum(0.f);
						}
						else if (std::isinf(alphaSum.y())) {
							LOG(ERROR) << StringPrintf(
								"Infinite alpha sum returned "
								"for pixel (%d, %d), sample %d. Setting to black.",
								pixel.x, pixel.y,
								(int)tileSampler->CurrentSampleNumber());
							alphaSum = Spectrum(0.f);
						}

                        const float mult = stats.hitCount > 0 ? float(stats.hitCount) / (stats.hitCount + stats.missCount) : 1.f;

						// Add the alpha sum directly to the image
                        filmTile->AddDirect(pixel, mult * alphaSum);
					}

					// Store extra data in exr channels
					if (moreDataLayers) {
						int offset = 0;
						for (int i = 0; i < techniques; ++i) {
							const Spectrum weightAvg = stats.weightsUpdates[i] > 0 ? stats.weights[i] / stats.weightsUpdates[i] : Spectrum(0.f);
							filmTile->AddData(pixel, weightAvg[0], offset++);
							filmTile->AddData(pixel, weightAvg[1], offset++);
							filmTile->AddData(pixel, weightAvg[2], offset++);
						}
						for (int i = 0; i < techniques; ++i) {
							const Spectrum alpha = alphasVector[0].GetAlpha(i, debug);
							filmTile->AddData(pixel, alpha[0], offset++);
							filmTile->AddData(pixel, alpha[1], offset++);
							filmTile->AddData(pixel, alpha[2], offset++);
						}
						const Spectrum mean = stats.contribsCount < 1 ? Spectrum(0.f) : stats.sumOfContribs / stats.contribsCount;
						const Spectrum var = stats.contribsCount < 2 ? Spectrum(0.f) :
							(stats.sumOfContribsSquared - stats.sumOfContribs * mean) / (stats.contribsCount - 1);
						const Spectrum varRel = (mean[0] == 0.f || mean[1] == 0.f || mean[2] == 0.f) ?
							Spectrum(0.f) : var / mean;
						filmTile->AddData(pixel, varRel[0], offset++);
						filmTile->AddData(pixel, varRel[1], offset++);
						filmTile->AddData(pixel, varRel[2], offset++);
					}
				}
				LOG(INFO) << "Finished image tile " << tileBounds;

				// Merge image tile into _Film_
                if (writeToImage) {
                    camera->film->MergeFilmTile(std::move(filmTile));
                }
				reporter.Update();
			}, nTiles);
			reporter.Done();
		}
		LOG(INFO) << writeToImage ? "Rendering finished" : "Dry run finished";

        const std::chrono::system_clock::time_point renderingEnd = std::chrono::system_clock::now();
        renderTime = std::chrono::duration_cast<std::chrono::milliseconds>(renderingEnd - renderingStart).count();
        if (writeToImage) {
            std::printf("\nRendering stats: samples %d, time %f s\n", sampler->samplesPerPixel, 0.001f * ((Float)renderTime));
            //statsFileTimes << 0.001f * ((Float)renderTime) << "\n";
            //statsFileSamples << sampler->samplesPerPixel << "\n";
        } else {
            std::printf("\nDry run stats: samples %d, time %f s\n", sampler->samplesPerPixel, 0.001f * ((Float)renderTime));
        }
	}

	// Save final image after rendering
    if (writeToImage) {
        camera->film->WriteImage();
    }

    //statsFileTimes.close();
    //statsFileSamples.close();
}

Spectrum OptimalMISIntegrator::Weight(const std::vector<Spectrum> &contribs, const std::vector<float> &pdfs, 
	const std::vector<int> &sampleCounts, AlphaEstimator& alphas, Stats *stats, bool debug) const {
	if (debug) {
		LOG(ERROR) << "Weighting";
		for (int i = 0; i < contribs.size(); ++i) {
			LOG(ERROR) << StringPrintf("contribs[%d]: ", i) << contribs[i].ToString() << "\n";
		}
		for (int i = 0; i < pdfs.size(); ++i) {
			LOG(ERROR) << StringPrintf("pdfs[%d]: %f\n", i, pdfs[i]);
		}
	}

	const int techniques = contribs.size();
	Spectrum L(0.f);
	if (techniques == 1) {
		if (!contribs[0].IsBlack()) {
			L += contribs[0] / pdfs[0];
		}
	}
	else {
		for (int i = 0; i < techniques; ++i) {
			if (misWeights == MISWeights::Optimal) {
				float denom = 0.f;
				Spectrum alphaComb(0.f);
				for (int j = 0; j < techniques; j++) {
					denom += sampleCounts[j] * pdfs[i * techniques + j];
					alphaComb += alphas.GetAlpha(j, debug) * pdfs[i * techniques + j];
				}
				if (pdfs[i * techniques + i] > 0.f) {
                    Spectrum wc = sampleCounts[i] / denom * (contribs[i] - alphaComb) + alphas.GetAlpha(i, debug);
					L += wc;

                    if (debug) {
                        LOG(ERROR) << StringPrintf("Weighted contrib: ") << wc.ToString() << "\n";
                    }

					if (contribs[i][0] != 0.f && contribs[i][1] != 0.f && contribs[i][2] != 0.f) {
                        for (int j = 0; j < techniques; j++) {
                            stats->weights[j] += sampleCounts[j] * pdfs[i * techniques + j] / denom * (Spectrum(1.f) - alphaComb / contribs[i]) + alphas.GetAlpha(j, debug) * pdfs[i * techniques + j] / contribs[i];
                            stats->weightsUpdates[j]++;
                        }
					}
				}
			}
			else if (!contribs[i].IsBlack()) {
				if (misWeights == MISWeights::Balance) {
					float denom = 0.f;
					for (int j = 0; j < techniques; j++) {
						denom += sampleCounts[j] * pdfs[i * techniques + j];
					}

					L += sampleCounts[i] / denom * contribs[i];

                    for (int j = 0; j < techniques; j++) {
                        stats->weights[j] += Spectrum(sampleCounts[j] * pdfs[i * techniques + j] / denom);
                        stats->weightsUpdates[j]++;
                    }
				}
				else { // misWeights == MISWeights::Power
					float denom = 0.f;
					for (int j = 0; j < techniques; j++) {
						denom += sampleCounts[j] * pdfs[i * techniques + j] * sampleCounts[j] * pdfs[i * techniques + j];
					}

					L += sampleCounts[i] / denom * contribs[i] * sampleCounts[i] * pdfs[i * techniques + i];

                    for (int j = 0; j < techniques; j++) {
                        stats->weights[j] += Spectrum(sampleCounts[j] * pdfs[i * techniques + j] * sampleCounts[j] * pdfs[i * techniques + j] / denom);
                        stats->weightsUpdates[j]++;
                    }
				}
			}
		}
	}

	stats->sumOfContribs += L;
	stats->sumOfContribsSquared += L * L;
	stats->contribsCount++;

    if (debug) {
        LOG(ERROR) << StringPrintf("Total weighted contrib: ") << L.ToString() << "\n";
    }

	return L;
}

}  // namespace pbrt
