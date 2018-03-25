#include "evoevo_def_state_pop.h"
#include "evoevo_def_base.h"
#include "evoevo_def_parameters.h"
#include "CFitnessEvaluator.h"
#include "CPrng.h"

#ifndef FITNESSRAND
#define FITNESSRAND
class CFitnessRandom : public CFitnessEvaluator{
	public:
		CFitnessRandom(CPrng* prng,TGeneElement maxNbClust,  TGeneElement maxNbDims);
		~CFitnessRandom();
		TFitness EvaluateFitness(TGenome* genomeToEvaluate, TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		TData* ClassifyData(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		TData* ComputePhenotype(TGenome* genomeToEvaluate, TIndexGeneElement geneSize);
};
#endif
