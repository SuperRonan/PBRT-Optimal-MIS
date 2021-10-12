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

		virtual void init(Scene const& scene) {};

		virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const = 0;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sammple) const = 0;

	};

	class LightSelector
	{
	protected:

		std::string strategy;
		std::shared_ptr<LightDistribution> distribution;
		std::unordered_map<const Light*, size_t> light_to_index;
		bool is_ns;

	public:
		
		std::vector<std::shared_ptr<Light>> lights;

		LightSelector(std::string const& strategy);

		void init(Scene const& scene);

		void init(Scene const& scene, std::vector<std::shared_ptr<Light>> const& lights);

		void select(const SurfaceInteraction& ref, Float lambda, const Light*& light, Float& pmf)const;

		Float pmf(const SurfaceInteraction& ref, const Light* light)const;
	};

	class GatheringTechnique : public LightSamplingTechnique
	{
	protected:

		static std::unordered_map<std::string, std::shared_ptr<LightSelector>> distributions;
		std::string strategy;
		std::shared_ptr<LightSelector> distribution;
		Scene const* scene;

	public:

		static void flushDistributions();

		std::shared_ptr<LightSelector> getSelector(std::vector<std::shared_ptr<Light>> const& lights, std::string const& specificity="");

		std::shared_ptr<LightSelector> getSelector(std::string const& specificity = "");

		static std::shared_ptr<LightSelector> getSelector(std::string const& strategy, Scene const& scene, std::vector<std::shared_ptr<Light>> const& lights, std::string const& specificity="");

		GatheringTechnique(std::string const& strategy);

		virtual void init(Scene const& scene) override;
	};

	class LiTechnique : public GatheringTechnique
	{
	public:

		LiTechnique(std::string const& strategy);

		virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const final override;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sample) const final override;
	};

	class GuidingTechnique : public GatheringTechnique
	{
	protected:

		GuidingDistribution::SamplingProjection projection_type;
		bool backup;

	public:

		GuidingTechnique(GuidingDistribution::SamplingProjection type, std::string const& strategy, bool backup=false);

		virtual void init(Scene const& scene) override;

		virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const final override;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sample) const final override;

	};


	//class LeTechnique : public GatheringTechnique
	//{
	//public:

	//	LeTechnique();

	//	virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const final override;

	//	virtual Float pdf(const SurfaceInteraction& ref, Sample const& sample) const final override;
	//};

	class SplattingTechnique : public LightSamplingTechnique
	{
	protected:
		Scene const* scene = nullptr;
		const Light* envmap = nullptr;
		float scene_radius = 0;
	public:
		SplattingTechnique();

		virtual void init(Scene const& scene) override;
	};

	class BSDFTechnique : public SplattingTechnique
	{
	public:

		BSDFTechnique();

		virtual void sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const final override;

		virtual Float pdf(const SurfaceInteraction& ref, Sample const& sample) const final override;
	};



}