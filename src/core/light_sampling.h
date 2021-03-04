#pragma once

#include <pbrt.h>
#include "light.h"
#include "lightdistrib.h"
#include "guiding.h"

namespace pbrt
{
	class LightSamplingTechnique
	{
	protected:

	public:

		enum class Type { Splatting, Gathering };

		const Type type;

		struct Sample
		{
			Type type;
			VisibilityTester vis;
			const Light* light = nullptr;
			Vector3f wi;
			Float pdf;
			Spectrum estimate;
			bool delta = false;
			bool visibility_passed = false;
		};

		LightSamplingTechnique(Type type);

		virtual void init(Scene const& scene, LightDistribution const& distrib) {};

		virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const = 0;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sammple) const = 0;

	};

	class GatheringTechnique : public LightSamplingTechnique
	{
	protected:

		LightDistribution const* distribution;
		std::unordered_map<const Light*, size_t> lightToIndex;
		Scene const* scene;

	public:

		GatheringTechnique();

		virtual void init(Scene const& scene, LightDistribution const& distrib) override;

		void selectLight(const SurfaceInteraction& ref, Float lambda, const Light*& light, Float& pdf) const;

		Float pdfSelectLight(const SurfaceInteraction& ref, const Light* light) const;
	};

	class LiTechnique : public GatheringTechnique
	{
	public:

		LiTechnique();

		virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const final override;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sample) const final override;
	};

	class LeTechnique : public GatheringTechnique
	{
	public:

		LeTechnique();

		virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const final override;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sample) const final override;
	};

	class SplattingTechnique : public LightSamplingTechnique
	{
	protected:
		Scene const* scene;
		const Light* envmap = nullptr;
		float scene_radius;
	public:
		SplattingTechnique();

		virtual void init(Scene const& scene, LightDistribution const&) override;
	};

	class BSDFTechnique : public SplattingTechnique
	{
	public:

		BSDFTechnique();

		virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const final override;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sample) const final override;
	};


	class GuidingTechnique : public GatheringTechnique
	{
	protected:

		std::vector<std::shared_ptr<Light>> lights;
		std::unique_ptr<LightDistribution> light_distrib;
		GuidingDistribution::SamplingProjection projection_type;

	public:

		GuidingTechnique(GuidingDistribution::SamplingProjection type);

		virtual void init(Scene const& scene, LightDistribution const&) override;

		virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const final override;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sample) const final override;

	};
}