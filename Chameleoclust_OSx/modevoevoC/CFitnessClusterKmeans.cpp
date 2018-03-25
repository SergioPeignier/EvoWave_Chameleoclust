#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "CFitnessClusterKmeans.h"

//creeer un parametre pour la taille max du cluster evaluator
CFitnessClusterKmeans::CFitnessClusterKmeans(CPrng* prng,TIndexGene maxNbGenes, TFitness nonCodingGenomeFitness, TNormExponent norm, TMutationLaw* statisticalLaw, TKmeansIterations nbMaxKmeansIterations,TGeneElement maxNbClust,  TGeneElement maxNbDims,  TGeneElement maxNbPoints):CFitnessEvaluator( maxNbClust,  maxNbDims){
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
	a_KMeanIterations                   = nbMaxKmeansIterations;
	a_modified_phenotype                = (TDataAggregator*)malloc(sizeof(TDataAggregator));
	a_phenotypeModifier                 = (TPoint*)malloc(sizeof(TPoint));
	a_phenotypeModifierlocal            = (TPoint*)malloc(sizeof(TPoint));
	AllocateDataAggregatorMem(a_modified_phenotype, maxNbClust, maxNbDims, maxNbPoints);
	AllocatePoint(a_phenotypeModifier,a_maxNbDims);
	AllocatePoint(a_phenotypeModifierlocal,a_maxNbDims);
}


void CFitnessClusterKmeans::AllocateDataAggregatorMem(TDataAggregator* dataAggregator, TGeneElement maxNbClust,  TGeneElement maxNbDims,  TGeneElement maxNbPoints){
	dataAggregator->p_arrayPoints = (TArrayPointsAggregator*) malloc(sizeof(TArrayPointsAggregator));
	dataAggregator->p_arrayPoints->p_elt = (TPointAggregator *) malloc(sizeof(TPointAggregator)*maxNbClust);
	for (TIndexPoint i = 0 ; i < maxNbClust; i++){
		dataAggregator->p_arrayPoints->p_elt[i].dimsArray  = (TDimAggregator * ) malloc(sizeof(TDimAggregator)*maxNbDims);
		for (TIndexPoint j = 0 ; j < maxNbDims; j++){
			dataAggregator->p_arrayPoints->p_elt[i].dimsArray[j].p_elt = (TCoordinateElement *) malloc(sizeof(TCoordinateElement)*maxNbPoints);
		}
	}
}
void CFitnessClusterKmeans::freeDataAggregator(TDataAggregator* dataAggregator){
	
	for (TIndexPoint i = 0 ; i < a_maxNbClust; i++){
		for (TIndexPoint j = 0 ; j < a_maxNbDims; j++){
			free(dataAggregator->p_arrayPoints->p_elt[i].dimsArray[j].p_elt);
		}
		free(dataAggregator->p_arrayPoints->p_elt[i].dimsArray);
	}
	free(dataAggregator->p_arrayPoints->p_elt);
	free(dataAggregator->p_arrayPoints);
	free(dataAggregator);
}

CFitnessClusterKmeans::~CFitnessClusterKmeans(){
	if (a_codingGenes->p_arrayGenes->nbElts>0) free(a_codingGenes->p_arrayGenes->p_elt);
	freeDataAggregator(a_modified_phenotype);
	free(a_codingGenes->p_arrayGenes);
	free(a_codingGenes);
	free(a_phenotypeModifier);
	free(a_phenotypeModifierlocal);
	}

//may be different for different norms
TFitness CFitnessClusterKmeans::DistanceNormalization( TFitness distance ,TGeneElement dimensionality,TNormExponent  norm){
	return distance*1.0/dimensionality;
	}


TFitness CFitnessClusterKmeans::AddtoDistance(TFitness additional, TNormExponent norm){
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


TFitness CFitnessClusterKmeans::AddtoDistanceNotInData(){

	switch (a_statisticalLaw->mutationLaw->law){
	
			case LAW_NORMAL:{

				return (TGeneElement)a_prng->gaussian(a_statisticalLaw->mutationLaw->mean, a_statisticalLaw->mutationLaw->standardDeviation);
				break;
				}
			case LAW_UNIFORM:{

				return a_prng->uniform(a_statisticalLaw->mutationLaw->min, a_statisticalLaw->mutationLaw->max);
				break;
				}

			}
			
	return 10000000;
	}

void CFitnessClusterKmeans::ForgetPhenotype(TData* phenotype){
	TIndexPoint i = 0;
	while (i < phenotype->p_arrayPoints->nbElts){
		phenotype->p_arrayPoints->p_elt[i].p_coordinates->nbElts = 0;
		phenotype->p_arrayPoints->p_elt[i].cluster = 0;
		phenotype->p_arrayPoints->p_elt[i].label   = 0;
		i++;
		}
	phenotype->p_arrayPoints->nbElts = 0;
	}


int CFitnessClusterKmeans::orderFunctionPhenotype(const void* A_i, const void* B_i){
	TGene * A = (TGene *) A_i;
	TGene * B = (TGene *) B_i;
  for(TGeneElement j = 0; j < INDEX_WEIG; j++){
    if(A->gene[j]>B->gene[j]) return 1;
	  if(A->gene[j]<B->gene[j]) return -1;
    }
	return 0;
	}

void CFitnessClusterKmeans::GetBoundaryGenesOfTypeX(TGenome* genomeSorted, TGeneElement type){
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

bool CFitnessClusterKmeans::HaveSameType_Cluster_Dimention(TGene gene1, TGene gene2){
	return gene1.gene[INDEX_TYPE] == gene2.gene[INDEX_TYPE] && gene1.gene[INDEX_CLUST] == gene2.gene[INDEX_CLUST] && gene1.gene[INDEX_DIM] == gene2.gene[INDEX_DIM];
	}

void  CFitnessClusterKmeans::CutAndSavePartOfGenome(TGenome* genomeSorted, TIndexGene start, TIndexGene stop){
	TIndexGene lengthToSave = stop - start;
	if (lengthToSave>0){
		memmove(&genomeSorted->p_arrayGenes->p_elt[0], &genomeSorted->p_arrayGenes->p_elt[start], lengthToSave*sizeof(TGene));
		genomeSorted->p_arrayGenes->nbElts = lengthToSave;
		}
	else{
		genomeSorted->p_arrayGenes->nbElts = 0;
		}
	}

void CFitnessClusterKmeans::TranscriptionGenome2CodingGenome(TGenome* genome, TGenome* genomeCoding){
	genomeCoding->p_arrayGenes->nbElts =  genome->p_arrayGenes->nbElts;
	memcpy(&genomeCoding->p_arrayGenes->p_elt[0], &genome->p_arrayGenes->p_elt[0], genome->p_arrayGenes->nbElts*sizeof(TGene));
	qsort(&genomeCoding->p_arrayGenes->p_elt[0], genome->p_arrayGenes->nbElts, sizeof(TGene) , orderFunctionPhenotype);
	GetBoundaryGenesOfTypeX(genomeCoding, CENTROID_GENE_TYPE);
	CutAndSavePartOfGenome(genomeCoding, a_start, a_stop);
	}

TGeneElement CFitnessClusterKmeans::AsociateGenesWeightedSum(TGenome* genomeCoding, TIndexGene start,TIndexGene  stop, TIndexGeneElement geneSize){
	TFitness position = 0;
	for (int i = start; i<=stop;i++){
		  position+=genomeCoding->p_arrayGenes->p_elt[i].gene[INDEX_POS];
		}
	return (TGeneElement)position;
	}


void CFitnessClusterKmeans::Coding2Clusters(TGenome* genomeCoding,TData* clusters,TIndexGeneElement geneSize){
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
				if (currentPoint > a_maxNbClust ) printf("POINT current %i  max %i \n", currentPoint,a_maxNbClust);
				if (currentCoord > a_maxNbDims )  printf("DIM current %i  max %i \n", currentCoord,a_maxNbDims);
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

void CFitnessClusterKmeans::InitializePhenotypeModified(){
	a_modified_phenotype->p_arrayPoints->nbElts = a_phenotype->p_arrayPoints->nbElts;
	for(TIndexPoint i = 0; i < a_modified_phenotype->p_arrayPoints->nbElts; i++){
		a_modified_phenotype->p_arrayPoints->p_elt[i].nbElts = a_phenotype->p_arrayPoints->p_elt[i].p_coordinates->nbElts;
		a_modified_phenotype->p_arrayPoints->p_elt[i].id = a_phenotype->p_arrayPoints->p_elt[i].label;
		for (TIndexPoint j = 0; j < a_modified_phenotype->p_arrayPoints->p_elt[i].nbElts; j ++){
			a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[j].nbElts = 0;
			a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[j].dimension = a_phenotype->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE];		
		}
	}
}


TData* CFitnessClusterKmeans::ComputePhenotype(TGenome* genomeToEvaluate,TIndexGeneElement geneSize){
	ForgetPhenotype(a_phenotype);
	TranscriptionGenome2CodingGenome(genomeToEvaluate,a_codingGenes);
	Coding2Clusters(a_codingGenes,a_phenotype,geneSize);
	InitializePhenotypeModified();
	return a_phenotype;
	}

TData* CFitnessClusterKmeans::ClassifyData(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	EvaluateFitness(genomeToEvaluate, geneSize, objectiveFunction, startObservationPosInData, stopObservationPosInData);
	return objectiveFunction;
	}

TFitness CFitnessClusterKmeans::AssignPoints(TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	if  (startObservationPosInData < stopObservationPosInData) return EvaluateFitnessInDataSubset(objectiveFunction, startObservationPosInData, stopObservationPosInData)*1.0/(stopObservationPosInData - startObservationPosInData);
	else return (EvaluateFitnessInDataSubset(objectiveFunction, 0, stopObservationPosInData) + EvaluateFitnessInDataSubset(objectiveFunction, startObservationPosInData, objectiveFunction->p_arrayPoints->nbElts))*1.0/(objectiveFunction->p_arrayPoints->nbElts);
}

TFitness CFitnessClusterKmeans::EvaluateFitness(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	ComputePhenotype(genomeToEvaluate,geneSize);
	if (a_phenotype->p_arrayPoints->nbElts == 0){
		return a_nonCodingGenomeFitness;
	}
	TKmeansIterations iterations = 0;
	TFitness result = 0;
	TBoolean convergence_k_means = FALSE;	
	while (!convergence_k_means && iterations < a_KMeanIterations){
		result = AssignPoints(objectiveFunction,startObservationPosInData,stopObservationPosInData);
		UpdatePhenotype();
		iterations ++;
	}
	return result;
}
    
TFitness CFitnessClusterKmeans::ComputeAllDistancesTo0(TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	TFitness  residualDistanceForCurrentCluster ;
	TFitness  sumMinDist = 0;
	for (TIndexPoint i = startObservationPosInData; i < stopObservationPosInData; i++){
		residualDistanceForCurrentCluster = 0;
		for (TIndexCoordinate j = 0; j < objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++){
			residualDistanceForCurrentCluster += AddtoDistance(objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE],a_norm);
			}
		sumMinDist += DistanceNormalization(residualDistanceForCurrentCluster,objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts,a_norm);
		}
	return sumMinDist*1.0;//(stopObservationPosInData-startObservationPosInData);
	}


//me imagino que aqui esta el problema
TFitness CFitnessClusterKmeans::EvaluateFitnessInDataSubset(TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData){
	
	TFitness sumMinDist = 0;

	TGeneElement     currentResidualDimensionality ;
	TIndexCoordinate currentClusterDim ;

	TFitness         minDistance2currentCluster;
	TFitness         minResidualDistanceWithCurrentCluster;

	TFitness         minLinealCombinationResidualExplainedCurrentCluster;
	TFitness         minLinealCombinationResidualExplained;

	TIndexCoordinate currentCoordinate ;

	
	TGeneElement     proximalCluster;

	for (TIndexPoint i = startObservationPosInData; i < stopObservationPosInData; i++){
		minLinealCombinationResidualExplained = FLT_MAX;
		proximalCluster                       = -1;
		//printf("_____________________\n");
		for (TIndexPoint j = 0; j < a_phenotype->p_arrayPoints->nbElts; j++){
			minResidualDistanceWithCurrentCluster               = 0;
			minDistance2currentCluster   		                = 0;
			currentCoordinate            		                = 0;
			currentClusterDim                                   = 0;
			minLinealCombinationResidualExplainedCurrentCluster = 0;
			currentResidualDimensionality                       = 0;
			a_phenotypeModifierlocal->p_coordinates->nbElts     = a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts;
			while(1){
				//Value in cluster and not in data , all data dims visited
				if (currentCoordinate == objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts && currentClusterDim < a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts){
					while (currentClusterDim < a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts){
						a_phenotypeModifierlocal->p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE] = a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE];
						a_phenotypeModifierlocal->p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE]  = a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE];
						minDistance2currentCluster += AddtoDistance(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE]-AddtoDistanceNotInData(),a_norm);
						currentClusterDim++;
					}
					break;
				}
				//Value in data and not in cluster , all cluster dims visited
				if (currentCoordinate < objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts && currentClusterDim == a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts){
					while (currentCoordinate < objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts){
						minResidualDistanceWithCurrentCluster += AddtoDistance(objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[POSITION_IN_COORDINATE],a_norm);
						currentResidualDimensionality++;
						currentCoordinate++;
					}
					break;
				}
				//all dims visited
				if (currentCoordinate == objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts && currentClusterDim == a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts){
					break;
				}
				//we have not finished visiting data dims not cluster dims
				if (currentCoordinate < objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->nbElts && currentClusterDim < a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->nbElts){
					//cluster and data share the same dim
					if(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE]   ==   objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[DIMENSION_IN_COORDINATE]){
						a_phenotypeModifierlocal->p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE] = a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE];
						a_phenotypeModifierlocal->p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE]  = objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[POSITION_IN_COORDINATE];
						minDistance2currentCluster += AddtoDistance(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE]-objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[POSITION_IN_COORDINATE],a_norm);
						currentCoordinate++;
						currentClusterDim++;
					}
					else{
						//cluster has a dim that data has not
						if(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE]   <    objectiveFunction->p_arrayPoints->p_elt[i].p_coordinates->p_elt[currentCoordinate].coordinate[DIMENSION_IN_COORDINATE]){
							a_phenotypeModifierlocal->p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE] = a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[DIMENSION_IN_COORDINATE];
							a_phenotypeModifierlocal->p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE]  = a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE];
							minDistance2currentCluster += AddtoDistance(a_phenotype->p_arrayPoints->p_elt[j].p_coordinates->p_elt[currentClusterDim].coordinate[POSITION_IN_COORDINATE]-AddtoDistanceNotInData(),a_norm);
							currentClusterDim++;
						}
						else{
							//data has a dim that cluster has not
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
			//printf("point %d --> cluster %d : distance: %lg\n",i,a_modified_phenotype->p_arrayPoints->p_elt[j].label,minLinealCombinationResidualExplainedCurrentCluster);
			if ( minLinealCombinationResidualExplainedCurrentCluster  < minLinealCombinationResidualExplained){
				/*
				printf("Modifier before being modified\n");
				for(TIndexCoordinate j = 0; j<a_phenotypeModifier->p_coordinates->nbElts; j++ ){
					printf("%d %d\n",a_phenotypeModifier->p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_phenotypeModifier->p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
				}
				printf("localModifier before modify\n");
				for(TIndexCoordinate j = 0; j<a_phenotypeModifierlocal->p_coordinates->nbElts; j++ ){
					printf("%d %d\n",a_phenotypeModifierlocal->p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_phenotypeModifierlocal->p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
				}
				*/
				ModifyPoint(a_phenotypeModifier,a_phenotypeModifierlocal);
				/*
				printf("Modifier after being modified\n");
				for(TIndexCoordinate j = 0; j<a_phenotypeModifier->p_coordinates->nbElts; j++ ){
					printf("%d %d\n",a_phenotypeModifier->p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_phenotypeModifier->p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
				}
				*/
				objectiveFunction->p_arrayPoints->p_elt[i].cluster	= a_phenotype->p_arrayPoints->p_elt[j].label;
				minLinealCombinationResidualExplained               = minLinealCombinationResidualExplainedCurrentCluster;
				proximalCluster = j;
    		}
		}
		// at this point we have visited all possible clusters for the given point
		//to comment
		//if (a_phenotypeModifier->p_coordinates->nbElts != a_modified_phenotype->p_arrayPoints->p_elt[proximalCluster].p_coordinates->nbElts) printf("BIG FAIL1");
		for (TIndexPoint k = 0; k < a_phenotypeModifier->p_coordinates->nbElts; k ++){
			a_modified_phenotype->p_arrayPoints->p_elt[proximalCluster].dimsArray[k].p_elt[a_modified_phenotype->p_arrayPoints->p_elt[proximalCluster].dimsArray[k].nbElts] = a_phenotypeModifier->p_coordinates->p_elt[k].coordinate[POSITION_IN_COORDINATE];			
			a_modified_phenotype->p_arrayPoints->p_elt[proximalCluster].dimsArray[k].nbElts ++;
			//to comment
			//if (a_modified_phenotype->p_arrayPoints->p_elt[proximalCluster].p_coordinates->p_elt[k].coordinate[DIMENSION_IN_COORDINATE] != a_phenotypeModifier->p_coordinates->p_elt[k].coordinate[DIMENSION_IN_COORDINATE]) printf("BIG FAIL2\n");
		}
		//printf("CHOOSE:::::point %d --> cluster %d : distance: %lg\n",i,a_modified_phenotype->p_arrayPoints->p_elt[proximalCluster].label,minLinealCombinationResidualExplained);
		sumMinDist +=minLinealCombinationResidualExplained;
	}
  return -sumMinDist;
}

void CFitnessClusterKmeans::UpdatePhenotype(){
	TBoolean convergence_k_means = TRUE;
	for (TIndexPoint i = 0; i < a_phenotype->p_arrayPoints->nbElts; i++){
		for (TIndexPoint k = 0; k < a_modified_phenotype->p_arrayPoints->p_elt[i].nbElts; k ++){
			if (a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].nbElts > 0){ 
				TCoordinateElement newValue = 0;
				if (a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].nbElts%2 == 1){
					newValue = qselect(a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].p_elt, a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].nbElts, a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].nbElts/2 + 1 );
				}
				if (a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].nbElts%2 == 0){
					TCoordinateElement valuerigth = qselect(a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].p_elt, a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].nbElts, a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].nbElts/2 );
					TCoordinateElement valueleft = qselect(a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].p_elt, a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].nbElts, a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].nbElts/2 + 1 );
					newValue = (valuerigth + valueleft )/2;
				}				
				convergence_k_means = convergence_k_means * (newValue == a_phenotype->p_arrayPoints->p_elt[i].p_coordinates->p_elt[k].coordinate[POSITION_IN_COORDINATE]);
				a_phenotype->p_arrayPoints->p_elt[i].p_coordinates->p_elt[k].coordinate[POSITION_IN_COORDINATE] = newValue;	
				a_modified_phenotype->p_arrayPoints->p_elt[i].dimsArray[k].nbElts = 0;
			}
		}
	}
}


void CFitnessClusterKmeans::PrintGenome(TGenome* genome){
	for(TIndexGene i = 0; i < genome->p_arrayGenes->nbElts; i++){
		for(TGeneElement j = 0; j < 5; j++){
			printf("%d,",genome->p_arrayGenes->p_elt[i].gene[j]);
			}
		printf("\n");
		}
	}

void CFitnessClusterKmeans::PrintPoint(TPoint point){
	for(TIndexCoordinate i = 0;i< point.p_coordinates->nbElts;i++){
		printf("(%d ,%d ); ",point.p_coordinates->p_elt[i].coordinate[DIMENSION_IN_COORDINATE], point.p_coordinates->p_elt[i].coordinate[POSITION_IN_COORDINATE]);
		}
	printf("\n");
	}

void  CFitnessClusterKmeans::PrintDataSet(TData* a_dataStream,TIndexPoint a_startObservationPosInData,TIndexPoint a_stopObservationPosInData){
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

void  CFitnessClusterKmeans::PrintDataSetToFile(FILE* f,TData* a_dataStream,TIndexPoint a_startObservationPosInData,TIndexPoint a_stopObservationPosInData){
	if  (a_startObservationPosInData < a_stopObservationPosInData){
		for (TIndexPoint i = a_startObservationPosInData; i<a_stopObservationPosInData;i++){
			fprintf(f,"point belonging to label: %d  and cluster: %d\n",a_dataStream->p_arrayPoints->p_elt[i].label, a_dataStream->p_arrayPoints->p_elt[i].cluster);
			for(TIndexCoordinate j = 0; j<a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				fprintf(f,"%d %d\n",a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
				}
			}
		}
	if (a_startObservationPosInData == a_stopObservationPosInData){
		for (TIndexPoint i = a_startObservationPosInData; i<a_dataStream->p_arrayPoints->nbElts;i++){
			fprintf(f,"point belonging to label: %d  and cluster: %d\n",a_dataStream->p_arrayPoints->p_elt[i].label, a_dataStream->p_arrayPoints->p_elt[i].cluster);
			for(TIndexCoordinate j = 0; j<a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				fprintf(f,"%d %d\n",a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
				}
			}

		for (TIndexPoint i = 0; i<a_startObservationPosInData;i++){
			fprintf(f,"point belonging to label: %d  and cluster: %d\n",a_dataStream->p_arrayPoints->p_elt[i].label, a_dataStream->p_arrayPoints->p_elt[i].cluster);
			for(TIndexCoordinate j = 0; j<a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				fprintf(f,"%d %d\n",a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
				}
			}
		}
	}
	
TCoordinateElement CFitnessClusterKmeans::qselect(TCoordinateElement *arr, TIndexPoint n, TIndexPoint k) {
  TIndexPoint i,ir,j,l,mid;
  TCoordinateElement a,temp;
  l=0;
  ir=n-1;
  for(;;) {
    if (ir <= l+1) { 
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir]);
      }
      return arr[k];
    }
    else {
      mid=(l+ir) >> 1; 
      SWAP(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
      }
      i=l+1; 
      j=ir;
      a=arr[l+1]; 
      for (;;) { 
	do i++; while (arr[i] < a); 
	do j--; while (arr[j] > a); 
	if (j < i) break; 
	SWAP(arr[i],arr[j]);
      } 
      arr[l+1]=arr[j]; 
      arr[j]=a;
      if (j >= k) ir=j-1; 
      if (j <= k) l=i;
    }
  }
}
