#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "CFitnessRandom.h"

CFitnessRandom::CFitnessRandom(CPrng* prng,TGeneElement maxNbClust,  TGeneElement maxNbDims):CFitnessEvaluator( maxNbClust,  maxNbDims){
	a_prng = prng;
	}


CFitnessRandom::~CFitnessRandom(){
	}
	
TFitness CFitnessRandom::EvaluateFitness(TGenome* genomeToEvaluate, TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	return a_prng->uniform();
	}
	
TData* CFitnessRandom::ClassifyData(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	return a_phenotype;
}


TData* CFitnessRandom::ComputePhenotype(TGenome* genomeToEvaluate, TIndexGeneElement geneSize){
	return a_phenotype;
}
