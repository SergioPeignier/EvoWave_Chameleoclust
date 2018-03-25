#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "CFitnessEvaluator.h"

CFitnessEvaluator::CFitnessEvaluator(TGeneElement maxNbClust, TGeneElement maxNbDims){
	a_maxNbClust =  maxNbClust;
	a_maxNbDims  =  maxNbDims;
	a_phenotype        = (TData*)malloc(sizeof(TData));
	AllocatePhenotypeMem(a_phenotype, a_maxNbClust, a_maxNbDims);
	a_prng             = NULL;
	}

CFitnessEvaluator::~CFitnessEvaluator(){
	freePhenotypeMem(a_phenotype);
}


void CFitnessEvaluator::freePhenotypeMem(TData* phenotype){
	if (phenotype != NULL){
		for (TIndexPoint i = 0;i<a_maxNbClust;i++){
			free(phenotype->p_arrayPoints->p_elt[i].p_coordinates->p_elt);
			free(phenotype->p_arrayPoints->p_elt[i].p_coordinates);
		}
		free(phenotype->p_arrayPoints->p_elt);
		free(phenotype->p_arrayPoints);
		free(phenotype);
	}
}

void CFitnessEvaluator::freePoint(TPoint* point){
	free(point->p_coordinates->p_elt);
	free(point->p_coordinates);
	free(point);
}

void CFitnessEvaluator::AllocatePoint(TPoint* point, TGeneElement max_array_length){
	point->p_coordinates = (TPointCoordinatesArray *)malloc(sizeof(TPointCoordinatesArray));
	point->p_coordinates->p_elt = (TCoordinate * )malloc(max_array_length*sizeof(TCoordinate));
	point->p_coordinates->nbElts=max_array_length;
	point->label = -1;
	point->cluster = -1;
	for (TIndexCoordinate j = 0; j<point->p_coordinates->nbElts; j++){
		point->p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE] = 0;
		point->p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]  = 0;
	}
}

void CFitnessEvaluator::ModifyPoint(TPoint* modified,TPoint* modifier){
	modified->label = modifier->label;
	modified->label = modifier->cluster;
	modified->p_coordinates->nbElts = modifier->p_coordinates->nbElts;
	for (TIndexCoordinate j = 0; j<modifier->p_coordinates->nbElts; j++){
		modified->p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE] = modifier->p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE];
		modified->p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]  = modifier->p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE];
	}
}

void CFitnessEvaluator::AllocatePhenotypeMem(TData* phenotype, TGeneElement maxNbClust, TGeneElement maxNbDims){
	phenotype->p_arrayPoints = (TArrayPoints*)malloc(sizeof(TArrayPoints));
	phenotype->p_arrayPoints->p_elt = (TPoint *)malloc(maxNbClust*sizeof(TPoint));
	phenotype->p_arrayPoints->nbElts = maxNbClust;
	for (TIndexPoint i = 0;i<maxNbClust;i++){
		phenotype->p_arrayPoints->p_elt[i].p_coordinates = (TPointCoordinatesArray *)malloc(sizeof(TPointCoordinatesArray));
		phenotype->p_arrayPoints->p_elt[i].p_coordinates->p_elt = (TCoordinate * )malloc(maxNbDims*sizeof(TCoordinate));
		phenotype->p_arrayPoints->p_elt[i].p_coordinates->nbElts=maxNbDims;
		phenotype->p_arrayPoints->p_elt[i].label = -1;
		phenotype->p_arrayPoints->p_elt[i].cluster = -1;
		for (TIndexCoordinate j = 0; j<maxNbDims; j++){
			phenotype->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE] = 0;
			phenotype->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]  = 0;
		}
	}
}

