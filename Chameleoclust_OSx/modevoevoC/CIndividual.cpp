#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <zlib.h>
#include <errno.h>
#include <string.h>
#include "CIndividual.h"

CIndividual::CIndividual(){
	a_fitness      = 0;
	a_fitnessIndex = 0;
    a_parentPosition = 0;
}

CIndividual::CIndividual(gzFile* backup_file , CPrng* prng){
	a_p_prng   =  prng;
	gzread((gzFile) backup_file, &a_position, sizeof(TIndexIndividual) );
	gzread((gzFile) backup_file, &a_fitness , sizeof(TFitness) );
	gzread((gzFile) backup_file, &a_fitnessIndex, sizeof(TIndexIndividual));
	gzread((gzFile) backup_file, &a_geneSize, sizeof(TIndexGeneElement));
	loadGenome(backup_file);
	}

void CIndividual::loadGenome(gzFile* backup_file){
	a_genome               =  (TGenome*) malloc(sizeof(TGenome));
	a_genome->p_arrayGenes =  (TArrayGenes*) malloc(sizeof(TArrayGenes));
	gzread((gzFile) backup_file, &a_genome->p_arrayGenes->maxNbElts, sizeof(TIndexGene));
	gzread((gzFile) backup_file, &a_genome->p_arrayGenes->nbElts, sizeof(TIndexGene));
	a_genome->p_arrayGenes->p_elt  =  (TGene *) malloc(a_genome->p_arrayGenes->maxNbElts*sizeof(TGene));
	gzread((gzFile) backup_file, &a_genome->p_arrayGenes->p_elt[0],a_genome->p_arrayGenes->nbElts*sizeof(TGene));
	}

void CIndividual::save( gzFile* backup_file ){
	gzwrite((gzFile) backup_file, &a_position, sizeof(TIndexIndividual) );
	gzwrite((gzFile) backup_file, &a_fitness , sizeof(TFitness) );
	gzwrite((gzFile) backup_file, &a_fitnessIndex, sizeof(TIndexIndividual));
	gzwrite((gzFile) backup_file, &a_geneSize, sizeof(TIndexGeneElement));
	saveGenome(backup_file);
	}

void CIndividual::saveGenome(gzFile* backup_file){
	gzwrite((gzFile) backup_file, &a_genome->p_arrayGenes->maxNbElts, sizeof(TIndexGene));
	gzwrite((gzFile) backup_file, &a_genome->p_arrayGenes->nbElts, sizeof(TIndexGene));
	gzwrite((gzFile) backup_file, &a_genome->p_arrayGenes->p_elt[0],a_genome->p_arrayGenes->nbElts*sizeof(TGene)); 
}

CIndividual::CIndividual(TIndexGene initNumGenes, TIndexGene maxNbGenes, TIndexGeneElement geneSize, TIndexIndividual position, CPrng* prng){
	a_p_prng        =  prng;
	a_position      =  position;
	a_fitness       =  0;
	a_fitnessIndex  =  0;
	a_genome        =  (TGenome*) malloc(sizeof(TGenome));
	a_geneSize      =  geneSize;
	a_genome->p_arrayGenes         =  (TArrayGenes*) malloc(sizeof(TArrayGenes));
	a_genome->p_arrayGenes->nbElts =  initNumGenes;
	a_genome->p_arrayGenes->maxNbElts = maxNbGenes;
	a_genome->p_arrayGenes->p_elt  =  (TGene *) malloc(maxNbGenes*sizeof(TGene));
	a_p_tempGen     = (TGene *) malloc(maxNbGenes * sizeof(TGene));
}

CIndividual::~CIndividual(){
	free(a_genome->p_arrayGenes->p_elt);
	free(a_genome->p_arrayGenes);
	free(a_genome);
	free(a_p_tempGen);
}

void CIndividual::CopyOtherIndividual(CIndividual* other){
	a_p_prng = other->a_p_prng;
	//a_position = other->a_position;
	a_fitness = other->a_fitness;
	a_parentPosition = other->a_position;
	memcpy(&a_genome->p_arrayGenes->p_elt[0],&other->a_genome->p_arrayGenes->p_elt[0],other->a_genome->p_arrayGenes->nbElts*sizeof(TGene));
	a_genome->p_arrayGenes->nbElts = other->a_genome->p_arrayGenes->nbElts;
}

void CIndividual::InitRandom(TMutationLaw** initialPositionsPDF, TBoundaryConditions** geneBoundaryConditions ){
	for(TIndexGene i = 0; i < a_genome->p_arrayGenes->nbElts; i++){
		for(TGeneElement j = 0; j < a_geneSize; j++){
			TGeneElement randomWalk = GeneElementRandomMutationWalk(initialPositionsPDF[j],0); //0 because we have to choose an initial position in the matrix, here we put 0, but should be changed
			a_genome->p_arrayGenes->p_elt[i].gene[j] = ComputeGeneElementAfterSubstitutionWithBoundary(0, randomWalk, geneBoundaryConditions[j]);
		}
	}
}

void CIndividual::BijectivePermutationGeneElements(TGene* pearl1, TGene* pearl2, TIndexGeneElement intergeniccutpoint, TGene* abstractGenesOrder){
	TGeneElement buffer = 0;
	for(TIndexGeneElement i = intergeniccutpoint; i< a_geneSize;i++){
		buffer = pearl1->gene[abstractGenesOrder->gene[i]];
		pearl1->gene[abstractGenesOrder->gene[i]] = pearl2->gene[abstractGenesOrder->gene[i]];
		pearl2->gene[abstractGenesOrder->gene[i]] = buffer;
	}
}
//spz, the end of 2 replaces the end of 1
void CIndividual::InjectivePermutationGeneElements(TGene* pearl1, TGene* pearl2, TIndexGeneElement intergeniccutpoint, TGene* abstractGenesOrder){
	for(TIndexGeneElement i = intergeniccutpoint; i< a_geneSize;i++){
		pearl1->gene[abstractGenesOrder->gene[i]] = pearl2->gene[abstractGenesOrder->gene[i]];
	}
}

TIndexGene CIndividual::ComputeTransposonSize(TIndexGene cutGeneStart, TIndexGene cutGeneEnd, TIndexGene len){
	if (cutGeneStart <= cutGeneEnd) return cutGeneEnd-cutGeneStart + 1;
	else return len - cutGeneStart + cutGeneEnd + 1;
}

void CIndividual::MakeTranslocationPermutation(TIndexGene cutGeneStart, TIndexGene cutGeneEnd,  TIndexGene insertionGene, TIndexGene len, TIndexGene transposonSize, TIntergenicCut* intergenicCutParameter){
	if ((len - transposonSize>=1)&&(intergenicCutParameter->cut)){
		TIndexGeneElement intergeniccutpoint = a_p_prng->uniform(intergenicCutParameter->mincutpoint, intergenicCutParameter->maxcutpoint);
		TIndexGene geneBeforeStart;
		geneBeforeStart = (cutGeneStart - 1) * (cutGeneStart > 0);
		BijectivePermutationGeneElements(&a_genome->p_arrayGenes->p_elt[geneBeforeStart],&a_genome->p_arrayGenes->p_elt[cutGeneEnd],intergeniccutpoint,intergenicCutParameter->genesAbstractOrder);
		BijectivePermutationGeneElements(&a_genome->p_arrayGenes->p_elt[insertionGene],&a_genome->p_arrayGenes->p_elt[cutGeneEnd],intergeniccutpoint,intergenicCutParameter->genesAbstractOrder);
	}
}

void CIndividual::MakeDeletionPermutation(TIndexGene cutGeneStart, TIndexGene cutGeneEnd, TIndexGene len, TIndexGene transposonSize, TIntergenicCut* intergenicCutParameter){
	if ((len - transposonSize>=0)&&(intergenicCutParameter->cut)){
		TIndexGeneElement intergeniccutpoint = a_p_prng->uniform(intergenicCutParameter->mincutpoint, intergenicCutParameter->maxcutpoint);
		TIndexGene geneBeforeStart;
		geneBeforeStart = (cutGeneStart - 1) * (cutGeneStart>0);
		BijectivePermutationGeneElements(&a_genome->p_arrayGenes->p_elt[geneBeforeStart],&a_genome->p_arrayGenes->p_elt[cutGeneEnd],intergeniccutpoint,intergenicCutParameter->genesAbstractOrder);
	}
}

void CIndividual::MakeDuplicationPermutation(TIndexGene cutGeneStart, TGene* cutGeneEndInTransposon, TIndexGene insertionGene,  TIndexGene len, TIndexGene transposonSize, TIntergenicCut* intergenicCutParameter){
	if ((len - transposonSize>=1)&&(intergenicCutParameter->cut)){
		TIndexGeneElement intergeniccutpoint = a_p_prng->uniform(intergenicCutParameter->mincutpoint, intergenicCutParameter->maxcutpoint);
		TIndexGene geneBeforeStart;
		geneBeforeStart = (cutGeneStart - 1) * (cutGeneStart>0);
		InjectivePermutationGeneElements(cutGeneEndInTransposon, &a_genome->p_arrayGenes->p_elt[geneBeforeStart], intergeniccutpoint, intergenicCutParameter->genesAbstractOrder);
		BijectivePermutationGeneElements(cutGeneEndInTransposon,&a_genome->p_arrayGenes->p_elt[insertionGene],intergeniccutpoint, intergenicCutParameter->genesAbstractOrder);
	}
}

void CIndividual::LargeDuplication(TIntergenicCut* intergenicCutParameter){
	if ( a_genome->p_arrayGenes->nbElts > 0 ){
		TIndexGene cutGeneStart=a_p_prng->uniform(0, a_genome->p_arrayGenes->nbElts);
		TIndexGene cutGeneEnd=a_p_prng->uniform(0, a_genome->p_arrayGenes->nbElts);
		TIndexGene insertionGene =a_p_prng->uniform(0, a_genome->p_arrayGenes->nbElts);
		LargeDuplication(cutGeneStart, cutGeneEnd, insertionGene, intergenicCutParameter);
	}
}

void CIndividual::LargeDuplication(TIndexGene cutGeneStart, TIndexGene cutGeneEnd, TIndexGene insertionGene, TIntergenicCut* intergenicCutParameter){
	TIndexGene len = a_genome->p_arrayGenes->nbElts;
	TIndexGene transposonSize = ComputeTransposonSize(cutGeneStart, cutGeneEnd, len);
	if (len + transposonSize < a_genome->p_arrayGenes->maxNbElts){
		if (cutGeneStart <= cutGeneEnd){
			memcpy(&a_p_tempGen[0], &a_genome->p_arrayGenes->p_elt[cutGeneStart] , transposonSize * sizeof(TGene));
		}
		else{
			memcpy(&a_p_tempGen[0],&a_genome->p_arrayGenes->p_elt[cutGeneStart] , (len - cutGeneStart) * sizeof(TGene) );
			memcpy(&a_p_tempGen[len - cutGeneStart] ,&a_genome->p_arrayGenes->p_elt[0] , (cutGeneEnd+1) * sizeof(TGene) );
		}
		MakeDuplicationPermutation(cutGeneStart, &a_p_tempGen[transposonSize - 1], insertionGene, len, transposonSize, intergenicCutParameter);		
		a_genome->p_arrayGenes->nbElts = a_genome->p_arrayGenes->nbElts + transposonSize;
		memmove( &a_genome->p_arrayGenes->p_elt[insertionGene + 1 + transposonSize], &a_genome->p_arrayGenes->p_elt[insertionGene + 1], (a_genome->p_arrayGenes->nbElts - insertionGene - transposonSize -1) * sizeof(TGene));
		memmove( &a_genome->p_arrayGenes->p_elt[insertionGene + 1], &a_p_tempGen[0], transposonSize * sizeof(TGene));
	}
}

TIndexGene CIndividual::LargeCut(TIndexGene cutGeneStart, TIndexGene cutGeneEnd){
	TIndexGene len = a_genome->p_arrayGenes->nbElts;
	TIndexGene newlength = 0;
	if (cutGeneStart <= cutGeneEnd){
		TIndexGene toCutLength = cutGeneEnd - cutGeneStart + 1;
		newlength = len - toCutLength;
		memmove(&a_genome->p_arrayGenes->p_elt[cutGeneStart], &a_genome->p_arrayGenes->p_elt[cutGeneEnd + 1], (len - cutGeneEnd - 1) * sizeof(TGene));
	}
	else{
		newlength = cutGeneStart - cutGeneEnd - 1;
		memmove(&a_genome->p_arrayGenes->p_elt[0], &a_genome->p_arrayGenes->p_elt[cutGeneEnd +1] , newlength  * sizeof(TGene));
	}
	return newlength;
}

void CIndividual::LargeDeletion(TIntergenicCut* intergenicCutParameter){
	if (a_genome->p_arrayGenes->nbElts > 0){
		TIndexGene cutGeneStart=a_p_prng->uniform(0, a_genome->p_arrayGenes->nbElts);
		TIndexGene cutGeneEnd=a_p_prng->uniform(0, a_genome->p_arrayGenes->nbElts);
		LargeDeletion(cutGeneStart, cutGeneEnd,intergenicCutParameter);
	}
}

void CIndividual::LargeDeletion(TIndexGene cutGeneStart, TIndexGene cutGeneEnd,TIntergenicCut* intergenicCutParameter){
	TIndexGene len = a_genome->p_arrayGenes->nbElts;
	TIndexGene transposonSize = ComputeTransposonSize(cutGeneStart, cutGeneEnd, len);
	MakeDeletionPermutation(cutGeneStart, cutGeneEnd, len, transposonSize, intergenicCutParameter);
	a_genome->p_arrayGenes->nbElts = LargeCut(cutGeneStart, cutGeneEnd);
}

void CIndividual::LargeTranslocation(TIntergenicCut* intergenicCutParameter){
	if (a_genome->p_arrayGenes->nbElts > 0){
		TIndexGene cutGeneStart = a_p_prng->uniform(0, a_genome->p_arrayGenes->nbElts);
		TIndexGene cutGeneEnd   = a_p_prng->uniform(0, a_genome->p_arrayGenes->nbElts);
		TIndexGene len          = a_genome->p_arrayGenes->nbElts;
		TIndexGene transposonSize = ComputeTransposonSize(cutGeneStart, cutGeneEnd, len);
		TIndexGene lengthAfterCut = len - transposonSize;
		if (lengthAfterCut > 0){
			TIndexGene insertionGene = a_p_prng->uniform(0, lengthAfterCut);
			LargeTranslocation(cutGeneStart, cutGeneEnd, insertionGene,intergenicCutParameter);
		}
	}
}

void CIndividual::LargeTranslocation(TIndexGene cutGeneStart, TIndexGene cutGeneEnd, TIndexGene insertionGene, TIntergenicCut* intergenicCutParameter){
	TIndexGene len = a_genome->p_arrayGenes->nbElts;
	TIndexGene transposonSize = ComputeTransposonSize(cutGeneStart, cutGeneEnd, len);
	MakeTranslocationPermutation(cutGeneStart,cutGeneEnd, insertionGene, len, transposonSize, intergenicCutParameter);
	if (cutGeneStart <= cutGeneEnd){
		memcpy(&a_p_tempGen[0], &a_genome->p_arrayGenes->p_elt[cutGeneStart] , transposonSize * sizeof(TGene));
		}
	else{
		memcpy(&a_p_tempGen[0],&a_genome->p_arrayGenes->p_elt[cutGeneStart] , (len - cutGeneStart) * sizeof(TGene) );
		memcpy(&a_p_tempGen[len - cutGeneStart] ,&a_genome->p_arrayGenes->p_elt[0] , (cutGeneEnd+1) * sizeof(TGene) );
	}
	LargeCut(cutGeneStart, cutGeneEnd);
	memmove( &a_genome->p_arrayGenes->p_elt[insertionGene + 1 + transposonSize], &a_genome->p_arrayGenes->p_elt[insertionGene + 1], (len - insertionGene - transposonSize -1) * sizeof(TGene));
	memmove( &a_genome->p_arrayGenes->p_elt[insertionGene + 1], &a_p_tempGen[0], transposonSize * sizeof(TGene));
}

void CIndividual::LargeInvertion(TIntergenicCut* intergenicCutParameter){
	if (a_genome->p_arrayGenes->nbElts > 0){
		TIndexGene cutGeneStart=a_p_prng->uniform(0, a_genome->p_arrayGenes->nbElts);
		TIndexGene cutGeneEnd=a_p_prng->uniform(0, a_genome->p_arrayGenes->nbElts);
		LargeInvertion(cutGeneStart, cutGeneEnd,intergenicCutParameter);
	}
}

void CIndividual::LargeInvertion(TIndexGene cutGeneStart, TIndexGene cutGeneEnd,TIntergenicCut* intergenicCutParameter){
	TGene temp;
	while (cutGeneStart > cutGeneEnd){
		temp = a_genome->p_arrayGenes->p_elt[cutGeneStart];
		a_genome->p_arrayGenes->p_elt[cutGeneStart] = a_genome->p_arrayGenes->p_elt[cutGeneEnd];
		a_genome->p_arrayGenes->p_elt[cutGeneEnd]   = temp;
		cutGeneStart++;
		cutGeneEnd--;
		if (cutGeneStart > (a_genome->p_arrayGenes->nbElts - 1)){
			cutGeneStart = 0;
		}
		if (cutGeneEnd < 0) {
			cutGeneEnd = a_genome->p_arrayGenes->nbElts - 1;
		}
	}
	while(cutGeneStart < cutGeneEnd){
		temp = a_genome->p_arrayGenes->p_elt[cutGeneStart];
		a_genome->p_arrayGenes->p_elt[cutGeneStart] = a_genome->p_arrayGenes->p_elt[cutGeneEnd];
		a_genome->p_arrayGenes->p_elt[cutGeneEnd] = temp;
		cutGeneStart++;
		cutGeneEnd--;
	}
}

void CIndividual::PrintGenome(){
	for(TIndexGene i = 0; i < a_genome->p_arrayGenes->nbElts; i++){
		for(TGeneElement j = 0; j < a_geneSize; j++){
			printf("%d,",a_genome->p_arrayGenes->p_elt[i].gene[j]);
			}
		printf("\n");
		}
	}
	
void CIndividual::Shuffle_genome(){
	if (a_genome->p_arrayGenes->nbElts>0) a_p_prng->shuffle((void*)a_genome->p_arrayGenes->p_elt,a_genome->p_arrayGenes->nbElts,sizeof(TGene));
	}

TGeneElement CIndividual::GeneElementRandomMutationWalk(TMutationLaw* mutationPDF, TIndexGeneElement index){
	switch (mutationPDF->mutationLaw->law){
			case LAW_NORMAL:{
				return (TGeneElement)a_p_prng->gaussian(mutationPDF->mutationLaw->mean, mutationPDF->mutationLaw->standardDeviation);
				break;
				}
			case LAW_UNIFORM:{
				return a_p_prng->uniform(mutationPDF->mutationLaw->min, mutationPDF->mutationLaw->max);
				break;
				}
			case LAW_TRANSITION_MATRIX:{
				return a_p_prng->wheelOfFotune(mutationPDF->mutationLaw->p_transitionMatrix->matrix[index],mutationPDF->mutationLaw->p_transitionMatrix->matrixSize);
				break;
				}
			}
	return 0;
	}

void CIndividual::DoGenePointSubstitution(TMutationLaw** mutationPDF,TBoundaryConditions** boundaryConditions ,TGeneEleMutationProbsLaw* geneEleMutationProbsLaw, TIndexGene mutationGene, TBoolean moreThanOnePointSubstitution){
	TIndexGeneElement obligatoryMutation = a_p_prng->wheelOfFotune(geneEleMutationProbsLaw->probGeneElementMutationWheelOfFortune, a_geneSize+1);
	a_genome->p_arrayGenes->p_elt[mutationGene].gene[obligatoryMutation] = ComputeGeneElementAfterSubstitutionWithBoundary(a_genome->p_arrayGenes->p_elt[mutationGene].gene[obligatoryMutation], GeneElementRandomMutationWalk(mutationPDF[obligatoryMutation], a_genome->p_arrayGenes->p_elt[mutationGene].gene[obligatoryMutation]), boundaryConditions[obligatoryMutation]);
	if (moreThanOnePointSubstitution){
		for(TIndexGeneElement i = 0; i < a_geneSize; i++){
			if (i!=obligatoryMutation){
				if (a_p_prng->uniform()< geneEleMutationProbsLaw->probGeneElementMutation[i] ){
					a_genome->p_arrayGenes->p_elt[mutationGene].gene[i] = ComputeGeneElementAfterSubstitutionWithBoundary(a_genome->p_arrayGenes->p_elt[mutationGene].gene[i], GeneElementRandomMutationWalk(mutationPDF[i], a_genome->p_arrayGenes->p_elt[mutationGene].gene[i]), boundaryConditions[i]);
				}
			}
		}
	}
}

TGeneElement CIndividual::ComputeGeneElementAfterSubstitutionWithBoundary(TGeneElement oldGeneElement, TGeneElement randomWalk, TBoundaryConditions* geneBoundaryConditions){
	TGeneElement newGeneElement = oldGeneElement + randomWalk;
	TGeneElement sizeBoundary = geneBoundaryConditions->max - geneBoundaryConditions->min - 1; //-1 = -2+1 : -2 excluding (bounds) +1 (to get the number of elements)
	if ((newGeneElement < geneBoundaryConditions->max)&&(newGeneElement > geneBoundaryConditions->min)) return(newGeneElement);
	return (sizeBoundary + (newGeneElement % sizeBoundary)) % sizeBoundary + geneBoundaryConditions->min + 1;
}

void CIndividual::DoAllLargeMutations(TSimulationParameters * simulationParameters){
	nbLargeMutations[LARGE_DUPLICATIONS_NB]   = a_p_prng->binomial( a_genome->p_arrayGenes->nbElts ,  simulationParameters->largeDuplicationRate);
	nbLargeMutations[LARGE_DELETIONS_NB]      = a_p_prng->binomial( a_genome->p_arrayGenes->nbElts ,  simulationParameters->largeDeletionRate);
	nbLargeMutations[LARGE_TRANSLOCATIONS_NB] = a_p_prng->binomial( a_genome->p_arrayGenes->nbElts ,  simulationParameters->largeTranslocationRate);
	nbLargeMutations[LARGE_INVERTIONS_NB]     = a_p_prng->binomial( a_genome->p_arrayGenes->nbElts ,  simulationParameters->largeInvertionRate);
	TIndexGene nbTotalLargeMutations = 0;
	TIndexGene random_value;
	for (TIndexGene i =0 ; i<NB_LARGE_MUTATIONS;i++){
		nbTotalLargeMutations += nbLargeMutations[i];
	}
	for ( TIndexGene i = nbTotalLargeMutations ; i > 0 ; i-- ){
		random_value = a_p_prng->uniform(0, i);
		if ( random_value < nbLargeMutations[LARGE_DUPLICATIONS_NB] ){
			LargeDuplication(simulationParameters->intergenicCutParameter);
			nbLargeMutations[LARGE_DUPLICATIONS_NB]--;
		}
		else if ( random_value < nbLargeMutations[LARGE_DUPLICATIONS_NB] + nbLargeMutations[LARGE_DELETIONS_NB] ) {
			LargeDeletion(simulationParameters->intergenicCutParameter);
			nbLargeMutations[LARGE_DELETIONS_NB]--;
		}
		else if ( random_value < nbLargeMutations[LARGE_DUPLICATIONS_NB] + nbLargeMutations[LARGE_DELETIONS_NB] + nbLargeMutations[LARGE_TRANSLOCATIONS_NB] ) {
			LargeTranslocation(simulationParameters->intergenicCutParameter);
			nbLargeMutations[LARGE_TRANSLOCATIONS_NB]--;
		}
		else {
			LargeInvertion(simulationParameters->intergenicCutParameter);
			nbLargeMutations[LARGE_INVERTIONS_NB]--;
		}
	}
}

void CIndividual::DoGenePointSubstitution(TMutationLaw** mutationPDF,TBoundaryConditions** boundaryConditions,TGeneEleMutationProbsLaw* geneEleMutationProbsLaw, TBoolean moreThanOnePointSubstitution){
	if (a_genome->p_arrayGenes->nbElts > 0){
		TIndexGene mutationGene = a_p_prng->uniform(0, a_genome->p_arrayGenes->nbElts);
		DoGenePointSubstitution( mutationPDF,boundaryConditions, geneEleMutationProbsLaw, mutationGene, moreThanOnePointSubstitution);
	}
}

void CIndividual::DoAllPointMutations(TSimulationParameters * simulationParameters){
	TIndexGene nbPointSubstitutions = a_p_prng->binomial( a_genome->p_arrayGenes->nbElts ,  simulationParameters->pointSubstitutionRate );
	for ( TIndexGene i = nbPointSubstitutions ; i > 0 ; i-- ){
		DoGenePointSubstitution(simulationParameters->mutationLaws,simulationParameters->boundaryConditions,simulationParameters->geneEleMutationProbsLaw, simulationParameters->moreThanOnePointSubstitution);
	}
}

void CIndividual::ComputeFitness(TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData ,CFitnessEvaluator * fitnessComputer ){
	a_fitness = fitnessComputer->EvaluateFitness(a_genome ,a_geneSize, dataStream, startObservationPosInData, stopObservationPosInData);
	}
	
TData* CIndividual::ComputePhenotype(CFitnessEvaluator * fitnessComputer){
	return fitnessComputer->ComputePhenotype(a_genome,a_geneSize);
}

void CIndividual::ComputeCodingRatio(){
	a_codingRatio = 0;
	if (a_genome->p_arrayGenes->nbElts != 0){
		for(TIndexGene i = 0; i < a_genome->p_arrayGenes->nbElts; i++){
			if (a_genome->p_arrayGenes->p_elt[i].gene[INDEX_TYPE] != PSEUDOGENE_GENE_TYPE){
				a_codingRatio++;
			}
		}
		a_codingRatio *= 1./a_genome->p_arrayGenes->nbElts;
	}
}

TIndividualFeatures CIndividual::getStatistics(){
	ComputeCodingRatio();
	TIndividualFeatures individualFeatures;
	individualFeatures.fitness         = a_fitness;
	individualFeatures.codingRatio     = a_codingRatio;
	individualFeatures.genomeSize      = a_genome->p_arrayGenes->nbElts;
	individualFeatures.parentIndex     = a_parentPosition;
	individualFeatures.individualIndex = a_position;
	return individualFeatures;
}

void CIndividual::SetGeneElement(TIndexGene geneIndex, TIndexGeneElement elementIndex, TGeneElement value){
	a_genome->p_arrayGenes->p_elt[geneIndex].gene[elementIndex] = value;
}

void CIndividual::SaveGenome(TFileName nameSaveFile){
	FILE *f;
    f = fopen(nameSaveFile, "w");
	for(TIndexGene i = 0; i < a_genome->p_arrayGenes->nbElts; i++){
		for(TGeneElement j = 0; j < a_geneSize; j++){
			fprintf(f,"%d,",a_genome->p_arrayGenes->p_elt[i].gene[j]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
}

void CIndividual::SaveFeatures(TFileName nameSaveFile, TGenerationsIndex curentGenerationIndex){
	TIndividualFeatures features = getStatistics();
	FILE *f;
    f = fopen(nameSaveFile, "a");
	fprintf(f,"%i %lg %lg %i %i %i\n",curentGenerationIndex,features.fitness,features.codingRatio,features.genomeSize,features.parentIndex,features.individualIndex);
	fclose(f);
}

void CIndividual::SavePhenotype( CFitnessEvaluator * fitnessComputer, TFileName nameSaveFile){
	TData* phenotype = ComputePhenotype(fitnessComputer);
	SaveDataSet(phenotype, 0, phenotype->p_arrayPoints->nbElts, nameSaveFile);
}

void CIndividual::SaveClassifiedData(TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData ,CFitnessEvaluator * fitnessComputer, TFileName nameSaveFile){
	TData* classifiedTable = fitnessComputer->ClassifyData(a_genome,a_geneSize, dataStream, startObservationPosInData, stopObservationPosInData);
	SaveDataSet(classifiedTable, startObservationPosInData, stopObservationPosInData, nameSaveFile);
}

void CIndividual::SaveDataSet(TData* dataset, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData, TFileName nameSaveFile){
	FILE *f;
    f = fopen(nameSaveFile, "w");
	if  (startObservationPosInData < stopObservationPosInData){
		for (TIndexPoint i = startObservationPosInData; i<stopObservationPosInData;i++){
			for(TIndexCoordinate j = 0; j<dataset->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				fprintf(f,"%d %d %d %d\n",dataset->p_arrayPoints->p_elt[i].cluster,dataset->p_arrayPoints->p_elt[i].label , dataset->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], dataset->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
			}
		}
	}
	if (startObservationPosInData == stopObservationPosInData){
		for (TIndexPoint i = 0; i<stopObservationPosInData;i++){
			for(TIndexCoordinate j = 0; j<=dataset->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				fprintf(f,"%d %d %d %d\n",dataset->p_arrayPoints->p_elt[i].cluster,dataset->p_arrayPoints->p_elt[i].label , dataset->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], dataset->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
			}
		}

		for (TIndexPoint i = startObservationPosInData; i<dataset->p_arrayPoints->nbElts;i++){
			for(TIndexCoordinate j = 0; j<=dataset->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				fprintf(f,"%d %d %d %d\n",dataset->p_arrayPoints->p_elt[i].cluster,dataset->p_arrayPoints->p_elt[i].label , dataset->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], dataset->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
			}
		}
	}
	fclose(f);
}

void CIndividual::ModifyPseudogenes(TMutationLaw** newPseudogenePDF, TBoundaryConditions** geneBoundaryConditions){
	for(TIndexGene i = 0; i < a_genome->p_arrayGenes->nbElts; i++){
		if (a_genome->p_arrayGenes->p_elt[i].gene[INDEX_TYPE] == 0){//replace with a more general operation
			for(TGeneElement j = 0; j < a_geneSize; j++){
				TGeneElement randomWalk = GeneElementRandomMutationWalk(newPseudogenePDF[j],0); //0 because we have to choose an initial position in the matrix, here we put 0, but should be changed
				a_genome->p_arrayGenes->p_elt[i].gene[j] = ComputeGeneElementAfterSubstitutionWithBoundary(0, randomWalk, geneBoundaryConditions[j]);
			}
		}
	}
}

int CIndividual::OrderFunctionGenotype(const void* A_i, const void* B_i){
	TGene * A = (TGene *) A_i;
	TGene * B = (TGene *) B_i;
	if(A->gene[INDEX_CLUST]>B->gene[INDEX_CLUST]) return 1;
	if(A->gene[INDEX_CLUST]<B->gene[INDEX_CLUST]) return -1;
	if(A->gene[INDEX_DIM]>B->gene[INDEX_DIM]) return 1;
	if(A->gene[INDEX_DIM]<B->gene[INDEX_DIM]) return -1;
	return 0;
	}

void CIndividual::OrderGenome(){
	qsort(&a_genome->p_arrayGenes->p_elt[0], a_genome->p_arrayGenes->nbElts, sizeof(TGene) , OrderFunctionGenotype);
}

