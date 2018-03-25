#include "evoevo_def_state_pop.h"
#include "evoevo_def_base.h"
#include "evoevo_def_parameters.h"
#include "CFitnessEvaluator.h"

#ifndef FITNESSTEST1
#define FITNESSTEST1
class CFitnessModeTest1 : public CFitnessEvaluator{
	public:
		TGeneElement* a_gapFunction;
		TGeneElement* a_activationPhenotype;
		TGeneElement* a_inhibitionPhenotype;
		TGeneElement* a_phenotype_array;
		TGeneElement a_max_array_length;
		TFitnessMode a_fitnessMode;
		TGeneElement width;
		TGeneElement top;
		TGeneElement height;
		TGeneElement weight;
		
		CFitnessModeTest1(TFitnessMode fitnessMode,TGeneElement maxNbClust,  TGeneElement maxNbDims);
		~CFitnessModeTest1();
		void ResetGapFunction( TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		TFitness EvaluateFitness(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		TGeneElement triangleLikeFunction(TGeneElement currentIndexInAxis,TGeneElement height, TGeneElement width , TGeneElement topTriangleIndex);
		void AddGeneElementsToArray(TGenome* genomeToEvaluate, TGeneElement* array);
		TData* ComputePhenotype(TGenome* genomeToEvaluate,TIndexGeneElement geneSize);
		TData* ClassifyData(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);

};
#endif
