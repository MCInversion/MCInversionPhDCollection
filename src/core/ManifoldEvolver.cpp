#include "ManifoldEvolver.h"

//
// ======================================================================================
//                    The strategy for 1D Curves in 2D space
// ---------------------------------------------------------------------------------------
//

void ManifoldCurveEvolutionStrategy::Preprocess()
{
}

void CustomManifoldCurveEvolutionStrategy::Preprocess()
{
}

void ManifoldCurveEvolutionStrategy::PerformEvolutionStep(unsigned int step)
{
}

bool ManifoldCurveEvolutionStrategy::ShouldRemesh()
{
	return false;
}

void ManifoldCurveEvolutionStrategy::Remesh()
{
}

void ManifoldCurveEvolutionStrategy::ExportCurrentState(unsigned int step)
{
}

void ManifoldCurveEvolutionStrategy::ExportFinalResult()
{
}

//
// ======================================================================================
//                    The strategy for 2D Surfaces in 3D space
// ---------------------------------------------------------------------------------------
//

void ManifoldSurfaceEvolutionStrategy::Preprocess()
{
}

void ManifoldSurfaceEvolutionStrategy::PerformEvolutionStep(unsigned int step)
{
}

bool ManifoldSurfaceEvolutionStrategy::ShouldRemesh()
{
	return false;
}

void ManifoldSurfaceEvolutionStrategy::Remesh()
{
}

void ManifoldSurfaceEvolutionStrategy::ExportCurrentState(unsigned int step)
{
}

void ManifoldSurfaceEvolutionStrategy::ExportFinalResult()
{
}

// ================================================
