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

	GatheringTechnique::GatheringTechnique() :
		LightSamplingTechnique(Type::Gathering)
	{}

	void GatheringTechnique::init(Scene const& scene, LightDistribution const& distrib)
	{
		distribution = &distrib;
		this->scene = &scene;
		for (size_t i = 0; i < scene.lights.size(); ++i)
		{
			lightToIndex[scene.lights[i].get()] = i;
		}
	}

	void GatheringTechnique::selectLight(const SurfaceInteraction& ref, Float lambda, const Light*& light, Float& pdf)const
	{
		const Distribution1D* distrib = distribution->Lookup(ref.p);
		int light_index = distrib->SampleDiscrete(lambda, &pdf);
		const std::shared_ptr<Light>& _light = scene->lights[light_index];
		light = _light.get();
	}

	Float GatheringTechnique::pdfSelectLight(const SurfaceInteraction& ref, const Light* light)const
	{
		const Distribution1D* distrib = distribution->Lookup(ref.p);
		const size_t light_id = lightToIndex.find(light)->second;
		Float light_pdf = distrib->DiscretePDF(light_id);
		return light_pdf;
	}

	LiTechnique::LiTechnique() :
		GatheringTechnique()
	{}

	void LiTechnique::sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const
	{
		Float light_pdf;
		selectLight(ref, lambda, sample.light, light_pdf);

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
			Float light_pdf = pdfSelectLight(ref, sample.light);
			// wi is not enough, what there are multiple points on the light in the direction of wi
			// It was enough for the previous use cases (no need to compute the PDF of zero contribution samples)
			Float p = sample.light->Pdf_Li(ref, sample.wi);
			return p * light_pdf;
		}
		return 0;
	}


	LeTechnique::LeTechnique() :
		GatheringTechnique()
	{}

	void LeTechnique::sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const
	{
		Float light_pdf;
		selectLight(ref, lambda, sample.light, light_pdf);
		Ray ray; Normal3f nlight; Float pdfDir, pdfA;
		sample.estimate = sample.light->Sample_Le(xi, xi, ref.time, &ray, &nlight, &pdfA, &pdfDir);
		// convert to solid angle density
		Float conversion = /*AbsDot(ray.d, nlight)*/ 1.0 / (ref.p - ray.o).LengthSquared();
		sample.pdf = pdfA * conversion;
		sample.pdf *= conversion;
		sample.wi = ray.d;
		Interaction it = Interaction(ray.o, ref.time, sample.light->mediumInterface);
		sample.vis = VisibilityTester(ref, it);
		sample.pdf *= light_pdf;
		if (!sample.estimate.IsBlack() && sample.pdf != 0)
		{
			sample.estimate *= ref.bsdf->f(ref.wo, sample.wi, BxDFType::BSDF_ALL) * AbsDot(sample.wi, ref.shading.n) / sample.pdf;
		}
		sample.type = this->type;
		sample.delta = sample.light->flags & (int(LightFlags::DeltaDirection) | int(LightFlags::DeltaPosition));
	}

	Float LeTechnique::pdf(const SurfaceInteraction& ref, Sample const& sample) const
	{
		// TODO
		Float light_pdf = pdfSelectLight(ref, sample.light);
		// wi is not enough, what there are multiple points on the light in the direction of wi
		// It was enough for the previous use cases (no need to compute the PDF of zero contribution samples)
		Float p = sample.light->Pdf_Li(ref, sample.wi);
		return p * light_pdf;
	}

	SplattingTechnique::SplattingTechnique() :
		LightSamplingTechnique(Type::Splatting)
	{}

	void SplattingTechnique::init(Scene const& scene, LightDistribution const&)
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



	GuidingTechnique::GuidingTechnique(GuidingDistribution::SamplingProjection type) :
		GatheringTechnique(),
		projection_type(type)
	{}

	void GuidingTechnique::init(Scene const& scene, LightDistribution const& distrib)
	{
		const auto acceptLight = [](Light const& light)
		{
			Point3f x, y, z;
			return light.GetTriangleVertices(&x, &y, &z);
		};

		distribution = &distrib;
		this->scene = &scene;

		for (size_t i = 0; i < scene.lights.size(); ++i)
		{
			const auto& light = scene.lights[i];
			if (acceptLight(*light))
			{
				lights.push_back(light);
				lightToIndex[light.get()] = lights.size() - 1;
			}
		}
		light_distrib = CreateLightSampleDistribution("power", scene, lights);
	}


	void GuidingTechnique::sample(const SurfaceInteraction& ref, Float lambda, Point2f const& xi, Sample& sample) const
	{
		sample.type = LightSamplingTechnique::Type::Gathering;
		if (lights.empty())
		{
			// Skip this sample, say that it failed
			sample.pdf = 0;
		}
		Float select_light_pmf;
		const Distribution1D* distrib = light_distrib->Lookup(ref.p);
		int light_index = distrib->SampleDiscrete(lambda, &select_light_pmf);
		sample.light = lights[light_index].get();

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
		else
		{
			// Skip this sample, say that it failed
			sample.pdf = 0;
		}
	}


	Float GuidingTechnique::pdf(const SurfaceInteraction& ref, Sample const& sample)const
	{
		if (lights.empty())	return 0;
		Float select_light_pmf = 0;
		const Distribution1D* distrib = light_distrib->Lookup(ref.p);
		if (lightToIndex.find(sample.light) == lightToIndex.end())	return 0;
		const int light_index = lightToIndex.find(sample.light)->second;
		select_light_pmf = distrib->DiscretePDF(light_index);

		GuidingDistribution gdstrb = GuidingDistribution(ref, *sample.light);
		if (!gdstrb.CanSample(projection_type))	return 0;

		Float pdf = gdstrb.Pdf(sample.wi, projection_type);
		return pdf * select_light_pmf;
	}
}

