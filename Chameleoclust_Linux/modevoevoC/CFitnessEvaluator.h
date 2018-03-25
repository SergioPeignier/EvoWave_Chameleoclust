#include "evoevo_def_state_pop.h"
#include "evoevo_def_base.h"
#include "evoevo_def_parameters.h"
#include "CPrng.h"
#ifndef FITNESSEVAL
#define FITNESSEVAL
class CFitnessEvaluator{
	public:
		TFitness a_fitness;
		TData* a_phenotype;
		TGeneElement a_maxNbClust;
		TGeneElement a_maxNbDims ;
		CPrng* a_prng;
		CFitnessEvaluator(TGeneElement maxNbClust, TGeneElement maxNbDims);
		virtual ~CFitnessEvaluator();
		virtual TFitness EvaluateFitness(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData) = 0;
		virtual TData* ComputePhenotype(TGenome* genomeToEvaluate, TIndexGeneElement geneSize) = 0;
		virtual void freePhenotypeMem(TData* phenotype);
		virtual void AllocatePhenotypeMem(TData* phenotype,TGeneElement maxNbClust, TGeneElement maxNbDims);
		virtual TData* ClassifyData(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData) = 0;
		virtual void freePoint(TPoint* point);
		virtual void AllocatePoint(TPoint* point, TGeneElement max_array_length);
		virtual void ModifyPoint(TPoint* modified,TPoint* modifier);
};
#endif
