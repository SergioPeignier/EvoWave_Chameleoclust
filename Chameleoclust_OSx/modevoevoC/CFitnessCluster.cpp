#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "CFitnessCluster.h"

CFitnessCluster::CFitnessCluster(CPrng* prng,TIndexGene maxNbGenes, TFitness nonCodingGenomeFitness, TNormExponent norm, TMutationLaw* statisticalLaw, TGeneElement maxNbClust,  TGeneElement maxNbDims):CFitnessEvaluator(  maxNbClust,   maxNbDims){
	a_maxNbClust = maxNbClust;
	a_maxNbDims = maxNbDims;
	a_prng    			             	= prng;
	a_nonCodingGenomeFitness        	= nonCodingGenomeFitness;
	a_codingGenes                       = (TGenome*) malloc(sizeof(TGenome));
	a_codingGenes->p_arrayGenes         = (TArrayGenes*) malloc(sizeof(TArrayGenes));
	a_codingGenes->p_arrayGenes->p_elt  = (TGene *) malloc(maxNbGenes*sizeof(TGene));
	a_codingGenes->p_arrayGenes->nbElts = 0;
	a_codingGenes->p_arrayGenes->maxNbElts = maxNbGenes;
	a_start                          	= 0;
	a_stop                           	= 0;
	a_norm                          	= norm;
	a_statisticalLaw                	= statisticalLaw;
}



CFitnessCluster::~CFitnessCluster(){
	free(a_codingGenes->p_arrayGenes->p_elt);
	free(a_codingGenes->p_arrayGenes);
	free(a_codingGenes);
}

TFitness CFitnessCluster::DistanceNormalization( TFitness distance ,TGeneElement dimensionality,TNormExponent  norm){
	return distance*1.0/dimensionality;
}


TFitness CFitnessCluster::AddtoDistance(TFitness additional, TNormExponent norm){
	if(norm == 1) return abs(additional);
	else {
    	TFitness ans  = 1;
    	TIndexPoint i = 0;
    	while (i < norm){
			ans *= additional;
      		i ++;
      	}
    return ans;
    }
}


TFitness CFitnessCluster::AddtoDistanceNotInData(){
	switch (a_statisticalLaw->mutationLaw->law){
		case LAW_NORMAL:{
			//printf("N %f %f\n",a_statisticalLaw->mutationLaw->mean, a_statisticalLaw->mutationLaw->standardDeviation);
			return (TGeneElement)a_prng->gaussian(a_statisticalLaw->mutationLaw->mean, a_statisticalLaw->mutationLaw->standardDeviation);
			break;
		}
		case LAW_UNIFORM:{
			//printf("U %f %f\n",a_statisticalLaw->mutationLaw->min*1.0, a_statisticalLaw->mutationLaw->max*1.0);
			return a_prng->uniform(a_statisticalLaw->mutationLaw->min, a_statisticalLaw->mutationLaw->max);
			break;
		}
	}
	return 1000;
}

void CFitnessCluster::ForgetPhenotype(){
	TIndexPoint i = 0;
	while (i < a_phenotype->p_arrayPoints->nbElts){
		a_phenotype->p_arrayPoints->p_elt[i].p_coordinates->nbElts = 0;
		i++;
	}
a_phenotype->p_arrayPoints->nbElts = 0;
}

int CFitnessCluster::orderFunctionPhenotype(const void* A_i, const void* B_i){
	TGene * A = (TGene *) A_i;
	TGene * B = (TGene *) B_i;
	for(TGeneElement j = 0; j < INDEX_WEIG; j++){	
		if(A->gene[j]>B->gene[j]) return 1;
		if(A->gene[j]<B->gene[j]) return -1;
    }
	return 0;
}

void CFitnessCluster::GetBoundaryGenesOfTypeX(TGenome* genomeSorted, TGeneElement type){
	TIndexGene i = 0;
	a_start = -1;
	a_stop  = -1;
	while (i < genomeSorted->p_arrayGenes->nbElts){
		if (a_start == -1 && genomeSorted->p_arrayGenes->p_elt[i].gene[INDEX_TYPE] == type) a_start = i;
		if (a_start != -1 && genomeSorted->p_arrayGenes->p_elt[i].gene[INDEX_TYPE] != type) {
			a_stop = i ;
			break;
		}
		i++;
	}
	if (a_start != -1 && a_stop == -1) a_stop= genomeSorted->p_arrayGenes->nbElts;
}

bool CFitnessCluster::HaveSameType_Cluster_Dimention(TGene gene1, TGene gene2){
	return gene1.gene[INDEX_TYPE] == gene2.gene[INDEX_TYPE] && gene1.gene[INDEX_CLUST] == gene2.gene[INDEX_CLUST] && gene1.gene[INDEX_DIM] == gene2.gene[INDEX_DIM];
}

void  CFitnessCluster::CutAndSavePartOfGenome(TGenome* genomeSorted, TIndexGene start, TIndexGene stop){
	TIndexGene lengthToSave = stop - start;
	if (lengthToSave>0){
		memmove(&genomeSorted->p_arrayGenes->p_elt[0], &genomeSorted->p_arrayGenes->p_elt[start], lengthToSave*sizeof(TGene));
		genomeSorted->p_arrayGenes->nbElts = lengthToSave;
	}
	else{
		genomeSorted->p_arrayGenes->nbElts = 0;
	}
}

void CFitnessCluster::TranscriptionGenome2CodingGenome(TGenome* genome, TGenome* genomeCoding){
	genomeCoding->p_arrayGenes->nbElts =  genome->p_arrayGenes->nbElts;
	memcpy(&genomeCoding->p_arrayGenes->p_elt[0], &genome->p_arrayGenes->p_elt[0], genome->p_arrayGenes->nbElts*sizeof(TGene));
	qsort(&genomeCoding->p_arrayGenes->p_elt[0], genome->p_arrayGenes->nbElts, sizeof(TGene) , orderFunctionPhenotype);
	GetBoundaryGenesOfTypeX(genomeCoding, CENTROID_GENE_TYPE);
	CutAndSavePartOfGenome(genomeCoding, a_start, a_stop);
}

/*
TGeneElement CFitnessCluster::AsociateGenesWeightedSum(TGenome* genomeCoding, TIndexGene start,TIndexGene  stop, TIndexGeneElement geneSize){
	TFitness position = 0;
	for (int i = start; i<=stop;i++){
    	position+=genomeCoding->p_arrayGenes->p_elt[i].gene[INDEX_POS];
	}
	return (TGeneElement)position;
}

TGeneElement CFitnessCluster::AsociateGenesWeightedSum(TGenome* genomeCoding, TIndexGene start,TIndexGene  stop, TIndexGeneElement geneSize){
	TFitness position = 0;
	for (int i = start; i<=stop;i++){
		if (abs(position) < abs(genomeCoding->p_arrayGenes->p_elt[i].gene[INDEX_POS])){
			position=genomeCoding->p_arrayGenes->p_elt[i].gene[INDEX_POS];
			}
	}
	return (TGeneElement)position;
}



TGeneElement CFitnessCluster::AsociateGenesWeightedSum(TGenome* genomeCoding, TIndexGene start,TIndexGene  stop, TIndexGeneElement geneSize){
	TFitness position = 0;
	for (int i = start; i<=stop;i++){
    	position+=genomeCoding->p_arrayGenes->p_elt[i].gene[INDEX_POS];
	}
	return (TGeneElement)position/(abs(stop-start)+0.00001);
}
*/

TGeneElement CFitnessCluster::AsociateGenesWeightedSum(TGenome* genomeCoding, TIndexGene start,TIndexGene  stop, TIndexGeneElement geneSize){
	TFitness position = 0;
	for (int i = start; i<=stop;i++){
    	position+=genomeCoding->p_arrayGenes->p_elt[i].gene[INDEX_POS];
	}
	return (TGeneElement)position;
}

void CFitnessCluster::Coding2Clusters(TGenome* genomeCoding,TData* clusters,TIndexGeneElement geneSize){
	TIndexGene actualEquals       = 0;
	TIndexPoint currentPoint      = 0;
	TIndexCoordinate currentCoord = 0;
	if (genomeCoding->p_arrayGenes->nbElts>0){
		TGeneElement currentCluster  = genomeCoding->p_arrayGenes->p_elt[0].gene[INDEX_CLUST];
		clusters->p_arrayPoints->p_elt[currentPoint].label = currentCluster;
		clusters->p_arrayPoints->nbElts++;
		for (TIndexGene i = 0; i < genomeCoding->p_arrayGenes->nbElts-1;i++){
			if (HaveSameType_Cluster_Dimention(genomeCoding->p_arrayGenes->p_elt[i],genomeCoding->p_arrayGenes->p_elt[i+1])){
				actualEquals ++;
			}
			else{
				clusters->p_arrayPoints->p_elt[currentPoint].p_coordinates->p_elt[currentCoord].coordinate[DIMENSION_IN_COORDINATE] =  genomeCoding->p_arrayGenes->p_elt[i].gene[INDEX_DIM];
				clusters->p_arrayPoints->p_elt[currentPoint].p_coordinates->p_elt[currentCoord].coordinate[POSITION_IN_COORDINATE]  =  AsociateGenesWeightedSum(genomeCoding, i - actualEquals, i, geneSize);
				currentCoord++;
				clusters->p_arrayPoints->p_elt[currentPoint].p_coordinates->nbElts++;
				if (currentCluster != genomeCoding->p_arrayGenes->p_elt[i+1].gene[INDEX_CLUST]){
					clusters->p_arrayPoints->nbElts++;
					currentCoord = 0;
					currentPoint++;
					currentCluster   = genomeCoding->p_arrayGenes->p_elt[i+1].gene[INDEX_CLUST];
					clusters->p_arrayPoints->p_elt[currentPoint].label = currentCluster;
				}
				actualEquals = 0;
			}
		}
		clusters->p_arrayPoints->p_elt[currentPoint].p_coordinates->p_elt[currentCoord].coordinate[DIMENSION_IN_COORDINATE] =  genomeCoding->p_arrayGenes->p_elt[genomeCoding->p_arrayGenes->nbElts-1].gene[INDEX_DIM];
		clusters->p_arrayPoints->p_elt[currentPoint].p_coordinates->p_elt[currentCoord].coordinate[POSITION_IN_COORDINATE]  =  AsociateGenesWeightedSum(genomeCoding, genomeCoding->p_arrayGenes->nbElts-1 - actualEquals, genomeCoding->p_arrayGenes->nbElts-1, geneSize);
		clusters->p_arrayPoints->p_elt[currentPoint].p_coordinates->nbElts++;
	}
}

TData* CFitnessCluster::ComputePhenotype(TGenome* genomeToEvaluate,TIndexGeneElement geneSize){
	ForgetPhenotype();
	TranscriptionGenome2CodingGenome(genomeToEvaluate,a_codingGenes);
	Coding2Clusters(a_codingGenes,a_phenotype,geneSize);
	return a_phenotype;
}
	
TData* CFitnessCluster::ClassifyData(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	EvaluateFitness(genomeToEvaluate, geneSize, objectiveFunction, startObservationPosInData, stopObservationPosInData);
	return objectiveFunction;
}

TFitness CFitnessCluster::EvaluateFitness(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	ComputePhenotype(genomeToEvaluate,geneSize);
	if (a_phenotype->p_arrayPoints->nbElts == 0){
		return a_nonCodingGenomeFitness;
	}
	if  (startObservationPosInData < stopObservationPosInData) return EvaluateFitnessInDataSubset(objectiveFunction, startObservationPosInData, stopObservationPosInData)*1.0/(stopObservationPosInData - startObservationPosInData);
	else return (EvaluateFitnessInDataSubset(objectiveFunction, 0, stopObservationPosInData) + EvaluateFitnessInDataSubset(objectiveFunction, startObservationPosInData, objectiveFunction->p_arrayPoints->nbElts))*1.0/(objectiveFunction->p_arrayPoints->nbElts);
}

TFitness CFitnessCluster::ComputeAllDistancesTo0(TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	TFitness  residualDistanceForCurrentCluster ;
	TFitness  sumMinDist = 0;
	for (TIndexPoint i = startObservationPosInData; i < stopObservationPosInData; i++){
		residualDistanceForCurrentCluster = 0;
		for (TIndexCoordinate j = 0; j < objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++){
			residualDistanceForCurrentCluster += AddtoDistance(objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE],a_norm);
		}
		sumMinDist += DistanceNormalization(residualDistanceForCurrentCluster,objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts,a_norm);
	}
	return sumMinDist*1.0/(stopObservationPosInData-startObservationPosInData);
}


TFitness CFitnessCluster::EvaluateFitnessInDataSubset(TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){

	TFitness sumMinDist = 0;
	TGeneElement     currentResidualDimensionality ;
	TIndexCoordinate currentClusterDim ;
	TFitness         minDistance2currentCluster;
	TFitness         minResidualDistanceWithCurrentCluster;
	TFitness         minLinealCombinationResidualExplainedCurrentCluster;
	TFitness         minLinealCombinationResidualExplained;
	TIndexCoordinate currentCoordinate ;

	for (TIndexPoint i = startObservationPosInData; i < stopObservationPosInData; i++){
		minLinealCombinationResidualExplained = FLT_MAX;
		for (TIndexPoint j = 0; j < a_phenotype->p_arrayPoints->nbElts; j++){
			minResidualDistanceWithCurrentCluster               = 0;
			minDistance2currentCluster   		                = 0;
			currentCoordinate            		                = 0;
			currentClusterDim                                   = 0;
			minLinealCombinationResidualExplainedCurrentCluster = 0;
			currentResidualDimensionality                       = 0;
			while(1){
				if (currentCoordinate == objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts && currentClusterDim < a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts){
					while (currentClusterDim < a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts){
						minDistance2currentCluster += AddtoDistance(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE]-AddtoDistanceNotInData(),a_norm);
						currentClusterDim++;
					}
					break;
				}
				if (currentCoordinate < objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts && currentClusterDim == a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts){
					while (currentCoordinate < objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts){
						minResidualDistanceWithCurrentCluster += AddtoDistance(objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[POSITION_IN_COORDINATE],a_norm);
						currentResidualDimensionality++;
						currentCoordinate++;
					}
					break;
				}
				if (currentCoordinate == objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts && currentClusterDim == a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts){
					break;
				}
				if (currentCoordinate < objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts && currentClusterDim < a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts){
					if(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE]   ==   objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[DIMENSION_IN_COORDINATE]){
						minDistance2currentCluster += AddtoDistance(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE]-objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[POSITION_IN_COORDINATE],a_norm);
						currentCoordinate++;
						currentClusterDim++;
					}

					else{
						if(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE]   <    objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[DIMENSION_IN_COORDINATE]){
							minDistance2currentCluster += AddtoDistance(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE]-AddtoDistanceNotInData(),a_norm);
							currentClusterDim++;
						}

						else{
							if(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE]   >   objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[DIMENSION_IN_COORDINATE]){
								minResidualDistanceWithCurrentCluster += AddtoDistance(objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[POSITION_IN_COORDINATE],a_norm);
								currentResidualDimensionality++;
								currentCoordinate++;
							}
						}
					}
				}
			}
			minLinealCombinationResidualExplainedCurrentCluster = DistanceNormalization(minDistance2currentCluster+minResidualDistanceWithCurrentCluster,currentClusterDim+currentResidualDimensionality,a_norm);
			if ( minLinealCombinationResidualExplainedCurrentCluster  < minLinealCombinationResidualExplained){
				objectiveFunction->p_arrayPoints->p_elt[i].cluster	= a_phenotype->p_arrayPoints->p_elt[j].label;
				minLinealCombinationResidualExplained               = minLinealCombinationResidualExplainedCurrentCluster;		
          	}
		}
		sumMinDist += minLinealCombinationResidualExplained;
	}
  return -sumMinDist;
}



void CFitnessCluster::PrintGenome(TGenome* genome){
	for(TIndexGene i = 0; i < genome->p_arrayGenes->nbElts; i++){
		for(TGeneElement j = 0; j < 5; j++){
			printf("%d,",genome->p_arrayGenes->p_elt[i].gene[j]);
		}
		printf("\n");
	}
}

void CFitnessCluster::PrintPoint(TPoint point){
	for(TIndexCoordinate i = 0;i< point.p_coordinates->nbElts;i++){
		printf("(%d ,%d ); ",point.p_coordinates->p_elt[i].coordinate[DIMENSION_IN_COORDINATE], point.p_coordinates->p_elt[i].coordinate[POSITION_IN_COORDINATE]);
	}
	printf("\n");
}

void  CFitnessCluster::PrintDataSet(TData* a_dataStream,TIndexPoint a_startObservationPosInData,TIndexPoint a_stopObservationPosInData){
	if  (a_startObservationPosInData < a_stopObservationPosInData){
		for (TIndexPoint i = a_startObservationPosInData; i<a_stopObservationPosInData;i++){
			printf("point belonging to label: %d  and cluster: %d\n",a_dataStream->p_arrayPoints->p_elt[i].label, a_dataStream->p_arrayPoints->p_elt[i].cluster);
			for(TIndexCoordinate j = 0; j<a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				printf("%d %d\n",a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
			}
		}
	}
	if (a_startObservationPosInData == a_stopObservationPosInData){
		for (TIndexPoint i = a_startObservationPosInData; i<a_dataStream->p_arrayPoints->nbElts;i++){
			printf("point belonging to label: %d  and cluster: %d\n",a_dataStream->p_arrayPoints->p_elt[i].label, a_dataStream->p_arrayPoints->p_elt[i].cluster);
			for(TIndexCoordinate j = 0; j<a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				printf("%d %d\n",a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
			}
		}

		for (TIndexPoint i = 0; i<a_startObservationPosInData;i++){
			printf("point belonging to label: %d  and cluster: %d\n",a_dataStream->p_arrayPoints->p_elt[i].label, a_dataStream->p_arrayPoints->p_elt[i].cluster);
			for(TIndexCoordinate j = 0; j<a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				printf("%d %d\n",a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
			}
		}
	}
}
