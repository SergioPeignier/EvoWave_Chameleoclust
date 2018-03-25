#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "CFitnessModeTest1.h"
#define D_MAX_WEIGHT_4_CODING 181
CFitnessModeTest1::CFitnessModeTest1( TFitnessMode fitnessMode, TGeneElement maxNbClust,  TGeneElement maxNbDims):CFitnessEvaluator( maxNbClust,   maxNbDims){
	a_max_array_length = maxNbDims;
	a_fitnessMode         = fitnessMode;
	if (a_max_array_length > 0){
		a_gapFunction         = (TGeneElement*) malloc(a_max_array_length * sizeof(TGeneElement));
		a_phenotype_array     = (TGeneElement*) malloc(a_max_array_length * sizeof(TGeneElement));
		a_activationPhenotype = (TGeneElement*) malloc(a_max_array_length * sizeof(TGeneElement));
		a_inhibitionPhenotype = (TGeneElement*) malloc(a_max_array_length * sizeof(TGeneElement));
		}
	width  = 0;
	top    = 0;
	height = 0;
	weight = 0;
	}



CFitnessModeTest1::~CFitnessModeTest1(){
	if (a_max_array_length > 0){
		free(a_gapFunction);
		free(a_activationPhenotype);
		free(a_inhibitionPhenotype);
		free(a_phenotype_array);
	}
}

void CFitnessModeTest1::ResetGapFunction( TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	for(TIndexGeneElement i = 0; i<a_max_array_length;i++){
		a_gapFunction[i] = 0;
		}
	if  (startObservationPosInData < stopObservationPosInData){
		for (TIndexPoint i = startObservationPosInData; i<stopObservationPosInData;i++){
			for(TIndexCoordinate j = 0; j<objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				if ((objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE] < a_max_array_length) && (objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE] >=0 )){
					a_gapFunction[objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE]] = objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE];
				}
			}
		}
	}
	if (startObservationPosInData == stopObservationPosInData){
		for (TIndexPoint i = 0; i<stopObservationPosInData;i++){
			for(TIndexCoordinate j = 0; j<objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				if ((objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE] < a_max_array_length) && (objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE] >=0 )){
					a_gapFunction[objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE]] = objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE];
					}
				}
			}

		for (TIndexPoint i = startObservationPosInData; i<objectiveFunction->p_arrayPoints->nbElts;i++){
			for(TIndexCoordinate j = 0; j<objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				if ((objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE] < a_max_array_length) && (objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE] >=0 )){
					a_gapFunction[objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE]] = objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE];
					}
				}
			}
		}
	}


TFitness CFitnessModeTest1::EvaluateFitness(TGenome* genomeToEvaluate , TIndexGeneElement geneSize,TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	ResetGapFunction(objectiveFunction, startObservationPosInData, stopObservationPosInData);
	AddGeneElementsToArray(genomeToEvaluate, a_gapFunction);
	TFitness g = 0;
	for(TIndexGeneElement i = 0; i<a_max_array_length;i++){
		g+=abs(a_gapFunction[i]);
		}
	return -g;
	}

TData* CFitnessModeTest1::ComputePhenotype(TGenome* genomeToEvaluate,TIndexGeneElement geneSize){
	for(TIndexGeneElement i = 0; i<a_max_array_length;i++){
		a_phenotype_array[i] = 0;
		}
	a_phenotype->p_arrayPoints->nbElts = 1;
	AddGeneElementsToArray(genomeToEvaluate, a_phenotype_array);
	for (TIndexGeneElement i = 0; i<a_max_array_length;i++){
		a_phenotype->p_arrayPoints->p_elt[0].p_coordinates->p_elt[i].coordinate[0] = i;
		a_phenotype->p_arrayPoints->p_elt[0].p_coordinates->p_elt[i].coordinate[1] = a_phenotype_array[i];
		}
	return a_phenotype;
	}
	
TData* CFitnessModeTest1::ClassifyData(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	return ComputePhenotype(genomeToEvaluate,geneSize);
}


void CFitnessModeTest1::AddGeneElementsToArray(TGenome* genomeToEvaluate, TGeneElement* array){
	switch(a_fitnessMode){
		case FITNESS_MODE_TEST_1:{
			for( TIndexGene i = 0; i<genomeToEvaluate->p_arrayGenes->nbElts; i++){
				if (genomeToEvaluate->p_arrayGenes->p_elt[i].gene[INDEX_TYPE] == CENTROID_GENE_TYPE){
					top    = genomeToEvaluate->p_arrayGenes->p_elt[i].gene[ INDEX_DIM ];
					height = genomeToEvaluate->p_arrayGenes->p_elt[i].gene[ INDEX_POS];
					if((top < a_max_array_length) && (top >= 0)){
						array[ top ] -= height;
						}
					}
				}
			break;
			}
		case FITNESS_MODE_TEST_2:{
			for( TIndexGene i = 0; i<genomeToEvaluate->p_arrayGenes->nbElts; i++){
				if (genomeToEvaluate->p_arrayGenes->p_elt[i].gene[INDEX_TYPE] == CENTROID_GENE_TYPE){
					width  = genomeToEvaluate->p_arrayGenes->p_elt[i].gene[INDEX_CLUST];
					top    = genomeToEvaluate->p_arrayGenes->p_elt[i].gene[ INDEX_DIM ];
					height = genomeToEvaluate->p_arrayGenes->p_elt[i].gene[ INDEX_POS];
					for(TGeneElement k = -width; k <= width; k++){
						if ((top+k >= 0) && (top+k)<a_max_array_length){
							array[top+k] -= triangleLikeFunction(top+k,height,width,top);
							//if (array[top+k] > a_max_array_length) array[top+k] = a_max_array_length;
							//if (array[top+k] < -a_max_array_length) array[top+k] = -a_max_array_length;
							}
						}
					}
				}
			break;
			}
		case FITNESS_MODE_TEST_3:{
			TGeneElement localPhenotype = 0;
			for(TIndexGeneElement i = 0; i<a_max_array_length;i++){
				a_inhibitionPhenotype[i] = 0;
				a_activationPhenotype[i] = 0;
				}

			for( TIndexGene i = 0; i<genomeToEvaluate->p_arrayGenes->nbElts; i++){
				if (genomeToEvaluate->p_arrayGenes->p_elt[i].gene[INDEX_TYPE] == CENTROID_GENE_TYPE){
					width  = genomeToEvaluate->p_arrayGenes->p_elt[i].gene[INDEX_CLUST];
					top    = genomeToEvaluate->p_arrayGenes->p_elt[i].gene[INDEX_DIM ];
					height = genomeToEvaluate->p_arrayGenes->p_elt[i].gene[INDEX_POS];
					weight = genomeToEvaluate->p_arrayGenes->p_elt[i].gene[INDEX_WEIG];
					if (weight<D_MAX_WEIGHT_4_CODING){
						float e = 1.0 - (weight*1.0/(a_max_array_length+1));
						for(TGeneElement k = -width; k <= width; k++){
							if ((top+k >= 0) && (top+k)<a_max_array_length){
								if (height<0) a_inhibitionPhenotype[top+k] += triangleLikeFunction(top+k,int(height*e),width,top);
								else a_activationPhenotype[top+k] +=triangleLikeFunction(top+k,int(height*e),width,top);
								}
							}
						}
					}
				}

			for(TIndexGeneElement i = 0; i<a_max_array_length;i++){
				if (a_inhibitionPhenotype[i]<(-a_max_array_length)) a_inhibitionPhenotype[i]=-a_max_array_length;
				if (a_activationPhenotype[i]>a_max_array_length)  a_activationPhenotype[i]=a_max_array_length;
				localPhenotype = a_inhibitionPhenotype[i] + a_activationPhenotype[i];
				if (localPhenotype<0) localPhenotype = 0;
				array[i]-=localPhenotype;
				}
			break;
			}
		}
	}






TGeneElement CFitnessModeTest1::triangleLikeFunction(TGeneElement currentIndexInAxis,TGeneElement height, TGeneElement width , TGeneElement topTriangleIndex){
	TGeneElement relativePosition = 0;
	if (currentIndexInAxis > topTriangleIndex ) relativePosition = topTriangleIndex + width - currentIndexInAxis;
	else relativePosition = currentIndexInAxis - (topTriangleIndex-width);
	return round(height*1./width * relativePosition);
	}
