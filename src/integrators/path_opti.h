#pragma once

#include <pbrt.h>
#include "light.h"
#include "integrator.h"
#include "lightdistrib.h"


#ifdef _MSC_VER
// Optimal MIS uses Eigen, which at some point defines an Infinity variable
// It creates a conflict with pbrt's Infinity macro
// Btw, why is Infinity defined as a macro when _MSC_VER is defined,
// and not as a static constexpr Float, like in the other case???
#pragma push_macro("Infinity")
#undef Infinity
#endif
#include <MIS/Estimators.h>
#ifdef _MSC_VER
#pragma pop_macro("Infinity")
#endif

namespace pbrt
{
	class LightSamplingTechnique
	{
	protected:



	public:

		enum class Type {Splatting, Gathering};

		const Type type;

		struct Sample
		{
			Type type;
			VisibilityTester vis;
			const Light* light;
			Vector3f wi;
			Float pdf;
			Spectrum estimate;
		};

		LightSamplingTechnique(Type type);

		virtual void init(Scene const& scene, LightDistribution const& distrib) {};

		virtual void sample(const SurfaceInteraction & ref, Float lambda, Point2f const& xi, Sample& sample) const = 0;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sammple) const = 0;

	};

	class GatheringTechnique : public LightSamplingTechnique
	{
	protected:
		
		LightDistribution const * distribution;
		std::unordered_map<const Light*, size_t> lightToIndex;
		Scene const * scene;

	public:

		GatheringTechnique();

		virtual void init(Scene const& scene, LightDistribution const& distrib) override;

		void selectLight(const SurfaceInteraction& ref, Float lambda, const Light*& light, Float& pdf) const;

		Float pdfSelectLight(const SurfaceInteraction& ref, const Light* light) const;
	};

	// TODO make an intermediate Type GatheringTechnique
	class LiTechnique : public GatheringTechnique
	{
	protected:


	public:

		LiTechnique();

		virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const final override;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sample) const final override;

	};



	// Though its behaviour is very close to fit in a SamplerIntegrator like the classic Path
	// We need to do some estimator management per pixel, and consider a whole set of samples for a pixel
	// which is not possible only in the Li function (without some dirty things).
	// We could make this class fit in the SamplerIntegrator by adding something like a 
	// pixelPreprocess and pixelPostprocess functions.
	class PathOptiIntegrator : public Integrator
	{
	public:
		
		struct Technique
		{
			std::shared_ptr<LightSamplingTechnique> technique;
			int n = 1;
		};
	
	protected:


		const int maxDepth;
		
		const std::string lightSampleStrategy;
		std::unique_ptr<LightDistribution> lightDistrib;

		using Estimator = MIS::Estimator<Spectrum, Float>;
		using EstimatorPtr = Estimator*;

		MIS::Heuristic heuristic;
		// One vector of vector of estimator for each thread
		//					-> One vector of estimator for each depth
		mutable std::vector<std::vector<EstimatorPtr>> _estimators_buffer;

		std::vector<Technique> techniques;

		std::shared_ptr<const Camera> camera;

		std::shared_ptr<Sampler> sampler;
		const Bounds2i pixelBounds;

	public:

		PathOptiIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
			std::shared_ptr<Sampler> sampler,
			const Bounds2i& pixelBounds,
			MIS::Heuristic h,
			std::vector<Technique> const& techniques,
			const std::string& lightSampleStrategy = "spatial");

		virtual ~PathOptiIntegrator();

		void Preprocess(const Scene& scene, Sampler& sampler);

		void directLighting(SurfaceInteraction const& isect, Scene const& scene, MemoryArena& arena, Sampler& sampler, Estimator& estimator, Float* wbuffer) const;

		// Returns the contribution of not computed with MIS, i.e. directly visibible lights (depth = 0)
		Spectrum TracePath(const RayDifferential& ray, const Scene& scene, Sampler& sampler, MemoryArena& arena, EstimatorPtr* estimators, Float* wbuffer, Spectrum beta=1, int depth = 0) const;

		virtual void Render(const Scene& scene) final;
	};

	PathOptiIntegrator* CreatePathOptiIntegrator(const ParamSet& params, 
		std::shared_ptr<Sampler> sampler, std::shared_ptr<const Camera> camera);

}