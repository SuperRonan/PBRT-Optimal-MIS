#include "integrators/path_opti.h"
#include "parallel.h"
#include "scene.h"
#include "camera.h"
#include "progressreporter.h"
#include "paramset.h"

namespace pbrt
{
	PathOptiIntegrator::PathOptiIntegrator(
		int maxDepth,
		std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler,
		const Bounds2i& pixelBounds,
		MIS::Heuristic h,
		std::vector<Technique> const& techniques,
		const std::string& lightSampleStrategy,
		bool conservative) :
		maxDepth(maxDepth),
		lightSampleStrategy(lightSampleStrategy),
		heuristic(h),
		camera(camera),
		sampler(sampler),
		techniques(techniques),
		pixelBounds(pixelBounds),
		conservative(conservative)
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
		const int threads = MaxThreadIndex();
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
					MIS::EstimatorCreateInfo<Float> eci;
					eci.heuristic = heuristic;
					eci.N = N;
					estimator = MIS::createEstimator<Spectrum, Float>(eci); 
					for (int i = 0; i < N; ++i)
					{
						estimator->setSampleForTechnique(i, techniques[i].n);
					}
				}
			}, threads);

		for (Technique& tech : techniques)
		{
			tech.technique->init(scene);
		}
	}

	void PathOptiIntegrator::Render(const Scene& scene)
	{
		if (techniques.empty())
		{
			Warning("Path Opti Integrator: Cannot render without any techniques!");
			return;
		}
		std::cout << "Preprocessing..." << std::endl;
		Preprocess(scene, *sampler);
		std::cout << "Conservative: " << conservative << std::endl;

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
					for (EstimatorPtr& estimator : estimators)
						estimator->loop();

					// Free _MemoryArena_ memory from computing image value
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

			Estimator& estimator = *estimators[depth]; // Get the estimator for the current depth

			// Draw the samples from multiple techniques to estimate direct lighting
			// And feed the samples to the estimator
			if (conservative)
				directLighting<true>(isect, scene, beta, arena, sampler, estimator, wbuffer);
			else
				directLighting<false>(isect, scene, beta, arena, sampler, estimator, wbuffer);


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

	template <bool CONSERVATIVE>
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
				//LOG(INFO) << "Drawing a sample with technique " << i << " (" << typeid(technique).name() << ")\n";
				technique.sample(it, lambda, xi, sample);
				if (sample.pdf == 0)
				{
					// Failed (that happens (rarely) in PBRT)
					// skip it (consider it as a singular zero
					continue;
				}
				const bool is_conservative = CONSERVATIVE;
				if constexpr (!is_conservative)
					if (sample.estimate.IsBlack())	continue; // skip the zero contribution sample

				Spectrum estimate = beta * sample.estimate / Float(ni);
				if (sample.type == LightSamplingTechnique::Type::Gathering)
				{
					//LOG(INFO) << "Visibility test\n";
					// visibility test (for gathering techniques)
					Spectrum visibility = sample.vis.Tr(scene, sampler);
					sample.visibility_passed = !visibility.IsBlack();
					estimate *= visibility;
					
					if constexpr (!is_conservative)
						if (visibility.IsBlack())	continue; // skip the zero contribution sample 
				}
				// Use the weights buffer to tmporarily store the PDFs
				Float sum = sample.pdf * Float(ni);
				// Effective PDFs actually
				Float*& pdfs = wbuffer;
				pdfs[i] = sum;
				//LOG(INFO) << "Computing PDFs\n";
				for (int l = 0; l < N; ++l)
				{
					if (l != i)
					{
						bool other_type_is_zero = sample.delta ||
							(sample.type == LightSamplingTechnique::Type::Gathering && !sample.visibility_passed) ||
							(sample.type == LightSamplingTechnique::Type::Splatting && sample.light == nullptr);
						if (sample.type != techniques[l].technique->type && (other_type_is_zero))
							pdfs[l] = 0;
						else
						{
							Float ql = techniques[l].technique->pdf(it, sample) * techniques[l].n;
							pdfs[l] = ql;
							sum += ql;
						}
					}
				}
				//LOG(INFO) << "Normalizing weights\n";
				for (int l = 0; l < N; ++l)
				{
					wbuffer[l] = pdfs[l] / sum;
				}
				//LOG(INFO) << "Adding estimate to estimator (" << typeid(estimator).name() << ") ...\n";
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
		else if (h_name == "progressive")
			h = MIS::Heuristic::Progressive;
		else
			Error(std::string(std::string("Heuristic ") + h_name + " is not recognized!").c_str());

		std::string lightStrategy = params.FindOneString("lightsamplestrategy",
			"power");

		// Sampling techniques
		std::vector<PathOptiIntegrator::Technique> techs;

		//int n_Le = params.FindOneInt("Le", 0);
		//if (n_Le)
		//{
		//	PathOptiIntegrator::Technique leTech;
		//	leTech.n = n_Le;
		//	leTech.technique = std::make_shared<LeTechnique>();
		//	techs.push_back(leTech);
		//}

		int n_bsdf = params.FindOneInt("BSDF", 0);
		if (n_bsdf)
		{
			PathOptiIntegrator::Technique bsdfTech;
			bsdfTech.n = n_bsdf;
			bsdfTech.technique = std::make_shared<BSDFTechnique>();
			techs.push_back(bsdfTech);
		}


		const std::vector<std::string> available_strats = { "uniform", "power", "spatial", "NSP", "NSS"};

		int n_SS = params.FindOneInt("SS", 0);
		if (n_SS)
		{
			PathOptiIntegrator::Technique ssTech;
			ssTech.n = n_SS;
			ssTech.technique = std::make_shared<GuidingTechnique>(GuidingDistribution::SamplingProjection::SphereSimple, lightStrategy);
			techs.push_back(ssTech);
		}

		for (std::string const& strat : available_strats)
		{
			std::string tech_str = strat + "-SS";
			int n = params.FindOneInt(tech_str, 0);
			if (n)
			{
				PathOptiIntegrator::Technique strat_ss;
				strat_ss.n = n;
				strat_ss.technique = std::make_shared<GuidingTechnique>(GuidingDistribution::SamplingProjection::SphereSimple, strat);
				techs.push_back(strat_ss);
			}
		}


		int n_SP = params.FindOneInt("SP", 0);
		if (n_SP)
		{
			PathOptiIntegrator::Technique spTech;
			spTech.n = n_SP;
			spTech.technique = std::make_shared<GuidingTechnique>(GuidingDistribution::SamplingProjection::SpherePrecise, lightStrategy);
			techs.push_back(spTech);
		}

		for (std::string const& strat : available_strats)
		{
			std::string tech_str = strat + "-SP";
			int n = params.FindOneInt(tech_str, 0);
			if (n)
			{
				PathOptiIntegrator::Technique strat_sp;
				strat_sp.n = n;
				strat_sp.technique = std::make_shared<GuidingTechnique>(GuidingDistribution::SamplingProjection::SpherePrecise, strat);
				techs.push_back(strat_sp);
			}
		}


		int n_PP = params.FindOneInt("PP", 0);
		if (n_PP)
		{
			PathOptiIntegrator::Technique ppTech;
			ppTech.n = n_PP;
			ppTech.technique = std::make_shared<GuidingTechnique>(GuidingDistribution::SamplingProjection::ParallelPlane, lightStrategy);
			techs.push_back(ppTech);
		}

		for (std::string const& strat : available_strats)
		{
			std::string tech_str = strat + "-PP";
			int n = params.FindOneInt(tech_str, 0);
			if (n)
			{
				PathOptiIntegrator::Technique strat_pp;
				strat_pp.n = n;
				strat_pp.technique = std::make_shared<GuidingTechnique>(GuidingDistribution::SamplingProjection::ParallelPlane, strat);
				techs.push_back(strat_pp);
			}
		}


		int n_Li = params.FindOneInt("Li", 0);
		if (n_Li)
		{
			PathOptiIntegrator::Technique liTech;
			liTech.n = n_Li;
			liTech.technique = std::make_shared<LiTechnique>(lightStrategy);
			techs.push_back(liTech);
		}

		for (std::string const& strat : available_strats)
		{
			std::string tech_str = strat + "-Li";
			int n = params.FindOneInt(tech_str, 0);
			if (n)
			{
				PathOptiIntegrator::Technique strat_li;
				strat_li.n = n;
				strat_li.technique = std::make_shared<LiTechnique>(strat);
				techs.push_back(strat_li);
			}
		}


		bool conservative = params.FindOneBool("conservative", true);

		return new PathOptiIntegrator(maxDepth, camera, sampler, pixelBounds, h, techs, lightStrategy, conservative); 
	}
}