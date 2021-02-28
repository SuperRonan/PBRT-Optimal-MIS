#include "integrators/path_opti.h"
#include "parallel.h"
#include "scene.h"
#include "camera.h"
#include "progressreporter.h"
#include "paramset.h"

namespace pbrt
{
	LightSamplingTechnique::LightSamplingTechnique(Type type):
		type(type)
	{}
	
	GatheringTechnique::GatheringTechnique():
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
		const Distribution1D * distrib = distribution->Lookup(ref.p);
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

	SplattingTechnique::SplattingTechnique():
		LightSamplingTechnique(Type::Splatting)
	{}

	BSDFTechnique::BSDFTechnique():
		SplattingTechnique()
	{}

	void BSDFTechnique::init(Scene const& scene, LightDistribution const&)
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


	PathOptiIntegrator::PathOptiIntegrator(
		int maxDepth,
		std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler,
		const Bounds2i& pixelBounds,
		MIS::Heuristic h,
		std::vector<Technique> const& techniques,
		const std::string& lightSampleStrategy) :
		maxDepth(maxDepth),
		lightSampleStrategy(lightSampleStrategy),
		heuristic(h),
		camera(camera),
		sampler(sampler),
		techniques(techniques),
		pixelBounds(pixelBounds)
	{}

	PathOptiIntegrator::~PathOptiIntegrator()
	{
		for (int tid = 0; tid < _estimators_buffer.size(); ++tid)
		{
			std::vector<EstimatorPtr> & estimators = _estimators_buffer[tid];
			for (EstimatorPtr & estimator : estimators)
				if (estimator)
					delete estimator;
		}
	}

	void PathOptiIntegrator::Preprocess(const Scene& scene, Sampler& sampler)
	{
		lightDistrib = CreateLightSampleDistribution(lightSampleStrategy, scene);
		const int threads = NumSystemCores();
		// numtechs
		int N = techniques.size();
		// Allocate all the estimators
		_estimators_buffer = std::vector<std::vector<Estimator*>>(threads, std::vector<Estimator*>());
		ParallelFor([&](int tid)
			{
				std::vector<Estimator*> & estimators = _estimators_buffer[tid];
				estimators.resize(maxDepth);
				for (int j = 0; j < estimators.size(); ++j)
				{
					Estimator*& estimator = estimators[j];
					estimator = MIS::createEstimator<Spectrum, Float>(heuristic, N); 
					for (int i = 0; i < N; ++i)
					{
						estimator->setSampleForTechnique(i, techniques[i].n);
					}
				}
			}, threads);

		for (Technique& tech : techniques)
		{
			tech.technique->init(scene, *lightDistrib);
		}
	}

	void PathOptiIntegrator::Render(const Scene& scene)
	{
		if (techniques.empty())
		{
			Warning("Path Opti Integrator: Cannot render without any techniques!");
			return;
		}
		Preprocess(scene, *sampler);

		// Compute number of tiles, _nTiles_, to use for parallel rendering
		Bounds2i sampleBounds = camera->film->GetSampleBounds();
		Vector2i sampleExtent = sampleBounds.Diagonal();
		const int tileSize = 16;
		Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
			(sampleExtent.y + tileSize - 1) / tileSize);

		ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");

		Film& film = *camera->film;

		ParallelFor2D([&](Point2i tile) {
			// Render section of image corresponding to _tile_

			// Allocate _MemoryArena_ for tile
			MemoryArena arena;

			// Get sampler instance for tile
			int seed = tile.y * nTiles.x + tile.x;
			std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

			// Compute sample bounds for tile
			int x0 = sampleBounds.pMin.x + tile.x * tileSize;
			int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
			int y0 = sampleBounds.pMin.y + tile.y * tileSize;
			int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
			Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
			LOG(INFO) << "Starting image tile " << tileBounds;

			std::unique_ptr<FilmTile> filmTile =
				camera->film->GetFilmTile(tileBounds);

			// Loop over pixels in tile to render them
			for (Point2i pixel : tileBounds) {
				{
					ProfilePhase pp(Prof::StartPixel);
					tileSampler->StartPixel(pixel);
				}

				// Do this check after the StartPixel() call; this keeps
				// the usage of RNG values from (most) Samplers that use
				// RNGs consistent, which improves reproducability /
				// debugging.
				if (!InsideExclusive(pixel, pixelBounds))
					continue;

				// Get the array of estimators for the pixel, and reset it
				std::vector<EstimatorPtr> & estimators = _estimators_buffer[ThreadIndex];
				for (EstimatorPtr& estimator : estimators)
					estimator->reset();

				Spectrum L = 0;

				do {
					// Initialize _CameraSample_ for current sample
					CameraSample cameraSample =
						tileSampler->GetCameraSample(pixel);

					// Generate camera ray for current sample
					RayDifferential ray;
					Float rayWeight =
						camera->GenerateRayDifferential(cameraSample, &ray);
					ray.ScaleDifferentials(
						1 / std::sqrt((Float)tileSampler->samplesPerPixel));

					// Draw path sample
					if (rayWeight > 0)
					{
						const int N = techniques.size();
						Float* wbuffer = (Float*)arena.Alloc(sizeof(Float) * N);
						L += TracePath(ray, scene, *tileSampler, arena, estimators.data(), wbuffer, rayWeight);
					}
					// Free _MemoryArena_ memory from computing image sample
					// value
					arena.Reset();
				} while (tileSampler->StartNextSample());
				L /= Float(sampler->samplesPerPixel);
				// Get the estimators estimations
				for (EstimatorPtr& estimator : estimators)
					L += estimator->solve(sampler->samplesPerPixel);
				// Add camera ray's contribution to the image tile
				filmTile->AddSample(Point2f(pixel.x + 0.5, pixel.y + 0.5), L);
			}
			LOG(INFO) << "Finished image tile " << tileBounds;
			camera->film->MergeFilmTile(std::move(filmTile));
			reporter.Update();
			}, nTiles);
		reporter.Done();
		camera->film->WriteImage();
	}

	Spectrum PathOptiIntegrator::TracePath(const RayDifferential& _ray, const Scene& scene, Sampler& sampler, MemoryArena& arena, EstimatorPtr* estimators, Float* wbuffer, Spectrum beta, int depth)const
	{
		RayDifferential ray = _ray;

		Spectrum res = 0;
		bool specularBounce = false;
		int bounce = 0;
		// We assume depth starts at zero
		for (depth; depth < maxDepth; ++depth)
		{
			SurfaceInteraction isect;
			bool foundIntersection = scene.Intersect(ray, &isect);

			
			if (depth == 0)
			{
				// Directly viewing the light
				// Special case: get the emission
				Spectrum direct;
				if (foundIntersection)
				{
					direct = isect.Le(-ray.d);
				}
				else
				{
					for (const auto& light : scene.infiniteLights)
						direct += light->Le(ray);
				}
				res = beta * direct;
			}
			
			if (!foundIntersection)	break;
			
			isect.ComputeScatteringFunctions(ray, arena, true);
			
			if (!isect.bsdf)
			{
				ray = isect.SpawnRay(ray.d);
				--depth;
				continue;
			}

			const Distribution1D* distrib = lightDistrib->Lookup(isect.p);

			Estimator& estimator = *estimators[depth]; // Get the estimator for the current depth

			// Draw the samples from multiple techniques to estimate direct lighting
			// And feed the samples to the estimator
			directLighting(isect, scene, beta, arena, sampler, estimator, wbuffer);

			// Sample the BSDF to continue the path
			Vector3f wo = -ray.d, wi;
			Float pdf;
			BxDFType flags;
			Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, BSDF_ALL, &flags);
			if (f.IsBlack() || pdf <= 0)	break;
			Float cos_wi = AbsDot(wi, isect.shading.n);
			beta *= f * cos_wi / pdf;
			specularBounce = (flags & BSDF_SPECULAR) != 0;
			ray = isect.SpawnRay(wi);
		}

		return res;
	}


	void PathOptiIntegrator::directLighting(SurfaceInteraction const& it, Scene const& scene, Spectrum const& beta, MemoryArena& arena, Sampler& sampler, Estimator& estimator, Float* wbuffer) const 
	{
		using Sample = LightSamplingTechnique::Sample;
		const int N = techniques.size();
		for (int i = 0; i < N; ++i)
		{
			const LightSamplingTechnique& technique = *techniques[i].technique;
			const int ni = techniques[i].n;
			for (int j = 0; j < ni; ++j)
			{
				Point2f xi = sampler.Get2D();
				Float lambda = sampler.Get1D();
				Sample sample;
				technique.sample(it, lambda, xi, sample);
				if (sample.pdf == 0)
				{
					// Failed (that happens (rarely) in PBRT)
					// skip it (consider it as a singular zero
					continue;
				}
				Spectrum estimate = beta * sample.estimate / Float(ni);
				if (sample.type == LightSamplingTechnique::Type::Gathering)
				{
					// visibility test (for gathering techniques)
					Spectrum visibility = sample.vis.Tr(scene, sampler);
					sample.visibility_passed = !visibility.IsBlack();
					estimate *= visibility;
				}
				// Use the weights buffer to tmporarily store the PDFs
				Float sum = sample.pdf;
				Float*& pdfs = wbuffer;
				pdfs[i] = sample.pdf;
				for (int l = 0; l < N; ++l)
				{
					if (l != i)
					{
						bool other_type_is_zero =
							(sample.type == LightSamplingTechnique::Type::Gathering && !sample.visibility_passed) ||
							(sample.type == LightSamplingTechnique::Type::Splatting && sample.light == nullptr);
						if (sample.type != techniques[l].technique->type && (sample.delta || other_type_is_zero))
							pdfs[l] = 0;
						else
						{
							Float pdf = techniques[l].technique->pdf(it, sample);
							pdfs[l] = pdf;
							sum += pdf;
						}
					}
				}
				for (int l = 0; l < N; ++l)
				{
					wbuffer[l] = pdfs[l] / sum;
				}

				estimator.addEstimate(estimate, wbuffer, i);
			}
		}
	}


	PathOptiIntegrator* CreatePathOptiIntegrator(const ParamSet& params, std::shared_ptr<Sampler> sampler, std::shared_ptr<const Camera> camera)
	{
		int maxDepth = params.FindOneInt("maxdepth", 5);
		
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

		std::string lightStrategy = params.FindOneString("lightsamplestrategy",
			"power");

		// Sampling techniques
		std::vector<PathOptiIntegrator::Technique> techs;

		int n_Li = params.FindOneInt("Li", 0);
		if (n_Li)
		{
			PathOptiIntegrator::Technique liTech;
			liTech.n = n_Li;
			liTech.technique = std::make_shared<LiTechnique>();
			techs.push_back(liTech);
		}

		int n_Le = params.FindOneInt("Le", 0);
		if (n_Le)
		{
			PathOptiIntegrator::Technique leTech;
			leTech.n = n_Le;
			leTech.technique = std::make_shared<LeTechnique>();
			techs.push_back(leTech);
		}

		int n_bsdf = params.FindOneInt("BSDF", 0);
		if (n_bsdf)
		{
			PathOptiIntegrator::Technique bsdfTech;
			bsdfTech.n = n_bsdf;
			bsdfTech.technique = std::make_shared<BSDFTechnique>();
			techs.push_back(bsdfTech);
		}

		return new PathOptiIntegrator(maxDepth, camera, sampler, pixelBounds, h, techs, lightStrategy); 
	}
}