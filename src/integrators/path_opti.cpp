#include "integrators/path_opti.h"
#include "parallel.h"
#include "scene.h"
#include "camera.h"
#include "progressreporter.h"

namespace pbrt
{
	PathOptiIntegrator::PathOptiIntegrator(
		int maxDepth,
		std::shared_ptr<const Camera> camera,
		std::shared_ptr<Sampler> sampler,
		const Bounds2i& pixelBounds,
		MIS::Heuristic h,
		const std::string& lightSampleStrategy) :
		maxDepth(maxDepth),
		lightSampleStrategy(lightSampleStrategy),
		heuristic(h),
		camera(camera),
		sampler(sampler),
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
		int threads = NumSystemCores();
		// numtechs
		int N = 2;
		// Allocate all the estimators
		_estimators_buffer = std::vector<std::vector<Estimator*>>(threads, std::vector<Estimator*>());
		ParallelFor([&](int tid)
			{
				std::vector<Estimator*> & estimators = _estimators_buffer[tid];
				for (int j = 0; j < estimators.size(); ++j)
				{
					Estimator*& estimator = estimators[j];
					estimator = MIS::createEstimator<Spectrum, Float>(heuristic, N);
				}
			}, threads);
	}

	void PathOptiIntegrator::Render(const Scene& scene)
	{
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
						const int N = estimators[0]->numTechs();
						Float* wbuffer = (Float*)arena.Alloc(sizeof(Float) * N);
						L = TracePath(ray, scene, *sampler, arena, estimators.data(), wbuffer, rayWeight);
					}
					// Free _MemoryArena_ memory from computing image sample
					// value
					arena.Reset();
				} while (tileSampler->StartNextSample());

				// Get the estimators estimations
				for (EstimatorPtr& estimator : estimators)
					L += estimator->solve(sampler->samplesPerPixel);
				// Add camera ray's contribution to image
				film.AddSplat(Point2f(pixel), L);
			}
			LOG(INFO) << "Finished image tile " << tileBounds;

			reporter.Update();
			}, nTiles);
		reporter.Done();
	}

	Spectrum PathOptiIntegrator::TracePath(const RayDifferential& _ray, const Scene& scene, Sampler& sampler, MemoryArena& arena, EstimatorPtr* estimators, Float* wbuffer, Spectrum beta, int depth)const
	{
		RayDifferential ray = _ray;

		Spectrum res = 0;
		bool specularBounce = false;
		int bounce = 0;
		// We assume depth starts at zero
		for (depth; depth <= maxDepth; ++depth)
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
			
			if (!foundIntersection || depth >= maxDepth)break;
			
			isect.ComputeScatteringFunctions(ray, arena, true);
			
			if (!isect.bsdf)
			{
				ray = isect.SpawnRay(ray.d);
				--depth;
				continue;
			}

			const Distribution1D* distrib = lightDistrib->Lookup(isect.p);

			if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) > 0)
			{
				
			}

		}

		return res;
	}

}