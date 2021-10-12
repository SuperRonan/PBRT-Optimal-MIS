#include "light_sampling.h"
#include "scene.h"
#include "interaction.h"
#include "bssrdf.h"

#define protected public
#define private public
#include "lights/diffuse.h"
#undef protected
#undef private


namespace pbrt
{
	LightSamplingTechnique::LightSamplingTechnique(Type type) :
		type(type)
	{}

	LightSelector::LightSelector(std::string const& strategy):
		strategy(strategy),
		is_ns(strategy == "NSP" || strategy == "NSS")
	{
		if (is_ns)
		{
			if (this->strategy == "NSS")
				this->strategy = "spatial";
			else if(this->strategy == "NSP")
				this->strategy = "power";
		}
	}

	void LightSelector::init(Scene const& scene)
	{
		init(scene, scene.lights);
	}

	void LightSelector::init(Scene const& scene, std::vector<std::shared_ptr<Light>> const& lights)
	{
		this->lights = lights;
		distribution = CreateLightSampleDistribution(strategy, scene);
		for (size_t i = 0; i < lights.size(); ++i)
		{
			light_to_index[lights[i].get()] = i;
		}
	}

	void LightSelector::select(const SurfaceInteraction& ref, Float lambda, const Light*& light, Float& pmf)const
	{
		const Distribution1D* distrib = distribution->Lookup(ref.p);
		int light_index = 0;
		if (is_ns && distrib->Count() > 1)
		{
			static thread_local std::vector<Float> contrib;
			contrib = std::vector<Float>(distrib->func);
			int max_id = 0;
			for (int i = 1; i < contrib.size(); ++i)
			{
				if (contrib[i] > contrib[max_id])
					max_id = i;
			}
			contrib[max_id] = 0;
			Distribution1D ns_distrib(contrib.data(), contrib.size());
			light_index = ns_distrib.SampleDiscrete(lambda, &pmf);
		}
		else
		{
			light_index = distrib->SampleDiscrete(lambda, &pmf);
		}
		const std::shared_ptr<Light>& _light = lights[light_index];
		light = _light.get();
	}

	Float LightSelector::pmf(const SurfaceInteraction& ref, const Light* light)const
	{
		const Distribution1D* distrib = distribution->Lookup(ref.p);
		Float pmf = 0;
		if (is_ns & distrib->Count() > 1)
		{
			static thread_local std::vector<Float> contrib;
			contrib = std::vector<Float>(distrib->func);
			int max_id = 0;
			for (int i = 1; i < contrib.size(); ++i)
			{
				if (contrib[i] > contrib[max_id])
					max_id = i;
			}
			contrib[max_id] = 0;
			Distribution1D ns_distrib(contrib.data(), contrib.size());
			
			const auto f = light_to_index.find(light);
			if (f != light_to_index.end())
			{
				pmf = ns_distrib.DiscretePDF(f->second);
			}
		}
		else
		{
			const auto f = light_to_index.find(light);
			if (f != light_to_index.end())
			{
				pmf = distrib->DiscretePDF(f->second);
			}
		}
		return pmf;
	}

	std::unordered_map<std::string, std::shared_ptr<LightSelector>> GatheringTechnique::distributions;

	GatheringTechnique::GatheringTechnique(std::string const& strategy) :
		LightSamplingTechnique(Type::Gathering),
		distribution(nullptr),
		strategy(strategy)
	{}

	void GatheringTechnique::flushDistributions()
	{
		distributions.clear();
	}

	std::shared_ptr<LightSelector> GatheringTechnique::getSelector(std::string const& strategy, Scene const& scene, std::vector<std::shared_ptr<Light>> const& lights, std::string const& specificity)
	{
		const auto f = distributions.find(strategy + specificity);
		if (f == distributions.end())
		{
			std::shared_ptr<LightSelector> selector = std::make_shared<LightSelector>(strategy);
			selector->init(scene, lights);
			distributions[strategy + specificity] = selector;
			return selector;
		}
		else
		{
			return f->second;
		}
	}

	std::shared_ptr<LightSelector> GatheringTechnique::getSelector(std::vector<std::shared_ptr<Light>> const& lights, std::string const& specificity)
	{
		return getSelector(strategy, *scene, lights, specificity);
	}

	std::shared_ptr<LightSelector> GatheringTechnique::getSelector(std::string const& specificity)
	{
		return getSelector(strategy, *scene, scene->lights, specificity);
	}

	void GatheringTechnique::init(Scene const& scene)
	{
		this->scene = &scene;
		distribution = getSelector();
	}

	LiTechnique::LiTechnique(std::string const& strategy) :
		GatheringTechnique(strategy)
	{}

	void LiTechnique::sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const
	{
		Float light_pdf;
		distribution->select(ref, lambda, sample.light, light_pdf);

		sample.estimate = sample.light->Sample_Li(ref, xi, &sample.wi, &sample.pdf, &sample.vis);
		sample.pdf *= light_pdf;
		if (!sample.estimate.IsBlack() && sample.pdf != 0)
		{
			sample.estimate *= ref.bsdf->f(ref.wo, sample.wi, BxDFType::BSDF_ALL) * AbsDot(sample.wi, ref.shading.n) / sample.pdf;
		}
		sample.type = this->type;
		sample.delta = sample.light->flags & (int(LightFlags::DeltaDirection) | int(LightFlags::DeltaPosition));
	}

	Float LiTechnique::pdf(const SurfaceInteraction& ref, Sample const& sample) const
	{
		if (sample.light)
		{
			Float light_pdf = distribution->pmf(ref, sample.light);
			// wi is not enough, what there are multiple points on the light in the direction of wi
			// It was enough for the previous use cases (no need to compute the PDF of zero contribution samples)
			Float p = sample.light->Pdf_Li(ref, sample.wi);
			return p * light_pdf;
		}
		return 0;
	}


	GuidingTechnique::GuidingTechnique(GuidingDistribution::SamplingProjection type, std::string const& strategy, bool backup) :
		GatheringTechnique(strategy),
		projection_type(type),
		backup(backup)
	{}

	void GuidingTechnique::init(Scene const& scene)
	{
		this->scene = &scene;
		std::vector<std::shared_ptr<Light>> samplable_lights;
		const auto acceptLight = [](Light const& light)
		{
			Point3f x, y, z;
			return light.GetTriangleVertices(&x, &y, &z);
		};
		for (size_t i = 0; i < scene.lights.size(); ++i)
		{
			const auto& light = scene.lights[i];
			if (acceptLight(*light))
			{
				samplable_lights.push_back(light);
			}
		}

		distribution = getSelector(samplable_lights, "OnlyTriangles");
	}


	void GuidingTechnique::sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const
	{
		sample.type = LightSamplingTechnique::Type::Gathering;
		if (distribution->lights.empty())
		{
			// Skip this sample, say that it failed
			sample.pdf = 0;
			return;
		}
		Float select_light_pmf;
		distribution->select(ref, lambda, sample.light, select_light_pmf);

		GuidingDistribution gdstrb = GuidingDistribution(ref, *sample.light);

		bool is_ok = gdstrb.CanSample(projection_type);
		if (is_ok)
		{
			sample.wi = gdstrb.Sample_wi(xi, projection_type, &sample.pdf);
			if (sample.pdf > 0)
			{
				sample.pdf *= select_light_pmf;
				sample.delta = false;
				const DiffuseAreaLight* _light = dynamic_cast<const DiffuseAreaLight*>(sample.light);
				assert(_light != nullptr);
				// Get the contribution
				{
					Ray ray = ref.SpawnRay(sample.wi);
					Float tHit;
					SurfaceInteraction light_inter;
					bool intersected = _light->shape->Intersect(ray, &tHit, &light_inter);
					if (!intersected)
					{
						sample.pdf = 0;
						return;
					}
					sample.estimate = ref.bsdf->f(ref.wo, sample.wi) * _light->L(light_inter, -ray.d) * AbsDot(ref.shading.n, sample.wi) / sample.pdf;
					sample.vis = VisibilityTester(ref, light_inter);
				}
			}
		}
		else if (backup)
		{
			sample.estimate = sample.light->Sample_Li(ref, xi, &sample.wi, &sample.pdf, &sample.vis);
			sample.pdf *= select_light_pmf;
			if (!sample.estimate.IsBlack() && sample.pdf != 0)
			{
				sample.estimate *= ref.bsdf->f(ref.wo, sample.wi, BxDFType::BSDF_ALL) * AbsDot(sample.wi, ref.shading.n) / sample.pdf;
			}
			sample.type = this->type;
			sample.delta = false;
		}
		else
		{
			// Skip this sample, say that it failed
			sample.pdf = 0;
		}
	}


	Float GuidingTechnique::pdf(const SurfaceInteraction& ref, Sample const& sample)const
	{
		if (distribution->lights.empty())	return 0;
		Float select_light_pmf = distribution->pmf(ref, sample.light);
		if (select_light_pmf == 0)	return 0;

		GuidingDistribution gdstrb = GuidingDistribution(ref, *sample.light);
		Float pdf = 0;
		if (gdstrb.CanSample(projection_type))
		{
			pdf = gdstrb.Pdf(sample.wi, projection_type);
		}
		else if (backup)
		{
			pdf = sample.light->Pdf_Li(ref, sample.wi);
		}
		return pdf * select_light_pmf;
	}



	//LeTechnique::LeTechnique() :
	//	GatheringTechnique()
	//{}

	//void LeTechnique::sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const
	//{
	//	Float light_pdf;
	//	selectLight(ref, lambda, sample.light, light_pdf);
	//	Ray ray; Normal3f nlight; Float pdfDir, pdfA;
	//	sample.estimate = sample.light->Sample_Le(xi, xi, ref.time, &ray, &nlight, &pdfA, &pdfDir);
	//	// convert to solid angle density
	//	Float conversion = /*AbsDot(ray.d, nlight)*/ 1.0 / (ref.p - ray.o).LengthSquared();
	//	sample.pdf = pdfA * conversion;
	//	sample.pdf *= conversion;
	//	sample.wi = ray.d;
	//	Interaction it = Interaction(ray.o, ref.time, sample.light->mediumInterface);
	//	sample.vis = VisibilityTester(ref, it);
	//	sample.pdf *= light_pdf;
	//	if (!sample.estimate.IsBlack() && sample.pdf != 0)
	//	{
	//		sample.estimate *= ref.bsdf->f(ref.wo, sample.wi, BxDFType::BSDF_ALL) * AbsDot(sample.wi, ref.shading.n) / sample.pdf;
	//	}
	//	sample.type = this->type;
	//	sample.delta = sample.light->flags & (int(LightFlags::DeltaDirection) | int(LightFlags::DeltaPosition));
	//}

	//Float LeTechnique::pdf(const SurfaceInteraction& ref, Sample const& sample) const
	//{
	//	// TODO
	//	Float light_pdf = pdfSelectLight(ref, sample.light);
	//	// wi is not enough, what there are multiple points on the light in the direction of wi
	//	// It was enough for the previous use cases (no need to compute the PDF of zero contribution samples)
	//	Float p = sample.light->Pdf_Li(ref, sample.wi);
	//	return p * light_pdf;
	//}

	SplattingTechnique::SplattingTechnique() :
		LightSamplingTechnique(Type::Splatting)
	{}

	void SplattingTechnique::init(Scene const& scene)
	{
		this->scene = &scene;
		// Find the envmap, if exists
		for (const auto& l : scene.infiniteLights)
		{
			if (l->flags & int(LightFlags::Infinite))
			{
				this->envmap = l.get();
				break;
			}
		}
		scene_radius = scene.WorldBound().Diagonal().Length() * 2.0;
	}

	BSDFTechnique::BSDFTechnique() :
		SplattingTechnique()
	{}

	void BSDFTechnique::sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const
	{
		sample.estimate = 0;
		sample.type = this->type;
		BxDFType delta;
		Spectrum fs = ref.bsdf->Sample_f(ref.wo, &sample.wi, xi, &sample.pdf, BSDF_ALL, &delta);
		if (sample.pdf > 0 || !fs.IsBlack())
		{
			fs *= AbsDot(sample.wi, ref.shading.n) / sample.pdf;
			sample.delta = delta & BxDFType::BSDF_SPECULAR;

			SurfaceInteraction lightIsect;
			Ray ray = ref.SpawnRay(sample.wi);
			bool foundIntersection = scene->Intersect(ray, &lightIsect);
			if (foundIntersection)
			{
				sample.estimate = fs * lightIsect.Le(-sample.wi);
				sample.light = lightIsect.primitive->GetAreaLight();
				sample.vis = VisibilityTester(ref, lightIsect);
			}
			else if (envmap)
			{
				sample.estimate = fs * envmap->Le(ray);
				sample.light = envmap;
				sample.vis = VisibilityTester(ref, Interaction(ray(scene_radius), -ray.d, ray.time, MediumInterface()));
			}
		}
	}

	Float BSDFTechnique::pdf(const SurfaceInteraction& ref, Sample const& sample)const
	{
		return ref.bsdf->Pdf(ref.wo, sample.wi);
	}

}

