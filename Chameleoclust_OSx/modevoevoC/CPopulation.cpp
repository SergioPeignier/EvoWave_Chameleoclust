#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <zlib.h>
#include <errno.h>
#include <string.h>
#include "CPopulation.h"
#include "CIndividual.h"

CPopulation::CPopulation(TPopulationSize populationSize, TSelectionStrength selectionPressure, CPrng* prng, TGenerationsIndex paramCurentGenerationIndex){
	a_p_prng                       =  prng;
	a_paramCurentGenerationindex   =  paramCurentGenerationIndex;
	a_p_reproductionProbs          =  (TReproductionProb*)malloc(populationSize*sizeof(TReproductionProb));
	a_populationSize               =  populationSize;
	a_populationPresent.nbElts     =  populationSize;
	a_selectionPressure            =  selectionPressure;
	a_populationPresent.p_elt      =  (CIndividual ** ) malloc(populationSize*sizeof(CIndividual*));
	perIndividualOffspring         =  (TIndexIndividual*) malloc(a_populationPresent.nbElts * sizeof(TIndexIndividual));
	a_availablePositions           =  (TIndexIndividual*) malloc(a_populationPresent.nbElts * sizeof(TIndexIndividual));
	a_occupiedPositions            =  (TIndexIndividual*) malloc(a_populationPresent.nbElts * sizeof(TIndexIndividual));
	GenerateReproductionProbas();
}

CPopulation::CPopulation(gzFile* backup_file , CPrng* prng){
	a_p_prng =  prng;
	gzread((gzFile) backup_file, &a_paramCurentGenerationindex, sizeof(TGenerationsIndex) );
	gzread((gzFile) backup_file, &a_populationSize, sizeof(TPopulationSize ));
	gzread((gzFile) backup_file, &a_selectionPressure, sizeof(TSelectionStrength));
	gzread((gzFile) backup_file, &a_nbAvailablePositions, sizeof(TIndexIndividual));
	gzread((gzFile) backup_file, &a_nbOccupiedPositions,sizeof(TIndexIndividual));

	a_p_reproductionProbs                     = (TReproductionProb*)malloc(a_populationSize*sizeof(TReproductionProb));
	perIndividualOffspring         =  (TIndexIndividual*) malloc(a_populationSize * sizeof(TIndexIndividual));
	a_availablePositions           =  (TIndexIndividual*) malloc(a_populationSize * sizeof(TIndexIndividual));
	a_occupiedPositions            =  (TIndexIndividual*) malloc(a_populationSize * sizeof(TIndexIndividual));
	GenerateReproductionProbas();
	loadIndividualsArray(backup_file,prng);
}

void CPopulation::loadIndividualsArray( gzFile* backup_file,CPrng* prng){
	a_populationPresent.p_elt      =  (CIndividual ** ) malloc(a_populationSize*sizeof(CIndividual*));
	a_populationPresent.nbElts     =  a_populationSize;
	for (TIndexIndividual i = 0;i<a_populationSize;i++){
		a_populationPresent.p_elt[i] = new CIndividual(backup_file,prng);
	}
}

void CPopulation::save( gzFile* backup_file ){
	gzwrite((gzFile) backup_file, &a_paramCurentGenerationindex, sizeof(TGenerationsIndex) );
	gzwrite((gzFile) backup_file, &a_populationSize, sizeof(TPopulationSize ));
	gzwrite((gzFile) backup_file, &a_selectionPressure, sizeof(TSelectionStrength));
	gzwrite((gzFile) backup_file, &a_nbAvailablePositions, sizeof(TIndexIndividual));
	gzwrite((gzFile) backup_file, &a_nbOccupiedPositions,sizeof(TIndexIndividual));
	saveIndividualsArray(backup_file);
	}

void CPopulation::saveIndividualsArray(gzFile* backup_file){
	for (TIndexIndividual i = 0;i<a_populationSize;i++){
		a_populationPresent.p_elt[i]->save(backup_file);
	}
}



CPopulation::~CPopulation(){
	for (TIndexIndividual i = 0; i < a_populationSize;i++){
		delete a_populationPresent.p_elt[i];
	}
	delete  a_bestOldIdividual;
	free(a_populationPresent.p_elt);
	free(a_p_reproductionProbs);
	free(perIndividualOffspring);
	free(a_availablePositions);
	free(a_occupiedPositions);
}

void CPopulation::GenerateRandomPopulation(TIndexGene initNumberGenes, TIndexGene maxNbGenes,TIndexGeneElement geneSize,TMutationLaw** initialPositionsPDF, TBoundaryConditions** geneBoundaryConditions ){
  a_bestOldIdividual             = new CIndividual(initNumberGenes, maxNbGenes, geneSize,-1, a_p_prng);
	for (TIndexIndividual i = 0;i<a_populationSize;i++){
		a_populationPresent.p_elt[i] = new CIndividual(initNumberGenes, maxNbGenes, geneSize,i, a_p_prng);
		a_populationPresent.p_elt[i]->InitRandom(initialPositionsPDF,  geneBoundaryConditions );
		}
	TIndexIndividual best_init_indiv = 	a_p_prng->uniform(0,a_populationSize);
	a_bestOldIdividual->CopyOtherIndividual(a_populationPresent.p_elt[best_init_indiv]);
	a_bestOldIdividual->a_position = a_populationPresent.p_elt[best_init_indiv]->a_position;
	}

void CPopulation::GenerateRandomPopulation(TIndexGene initNumberGenes,TIndexGene maxNbGenes, TIndexGeneElement geneSize,TMutationLaw** initialPositionsPDF, TBoundaryConditions** geneBoundaryConditions , TIndexIndividual initGene, TIndexIndividual endGene){
  a_bestOldIdividual             = new CIndividual(initNumberGenes,maxNbGenes, geneSize,-1, a_p_prng);
  for (TIndexIndividual i = initGene;i<endGene;i++){
		a_populationPresent.p_elt[i] = new CIndividual(initNumberGenes,maxNbGenes, geneSize,i, a_p_prng);
		a_populationPresent.p_elt[i]->InitRandom(initialPositionsPDF , geneBoundaryConditions);
		}
	TIndexIndividual best_init_indiv = 	a_p_prng->uniform(0,a_populationSize);
	a_bestOldIdividual->CopyOtherIndividual(a_populationPresent.p_elt[best_init_indiv]);
	a_bestOldIdividual->a_position = a_populationPresent.p_elt[best_init_indiv]->a_position;
	}

void CPopulation::PrintAllGenomes(){

	printf("***--------------Population at generation %d -----------***\n",a_paramCurentGenerationindex);
	for (TIndexIndividual i = 0;i<a_populationSize;i++){
		printf("--------------Individual %d-----------\n",a_populationPresent.p_elt[i]->a_position);
		a_populationPresent.p_elt[i]->PrintGenome();
		}
	}

void CPopulation::SetPRNG(CPrng* new_prng){
	for (TIndexIndividual i = 0;i<a_populationSize;i++){
		a_populationPresent.p_elt[i]->a_p_prng = new_prng;
		}
	}

TMutationRate CPopulation::ReproductionProba(TIndexIndividual rank){
	//formula:     (c-1)*1.0/(c^N-1)*c^(N-r)
	return (a_selectionPressure-1)*1.0/(pow(a_selectionPressure,1.0*a_populationSize)-1)*pow(a_selectionPressure,1.*(a_populationSize-rank));
	}

void CPopulation::GenerateReproductionProbas(){
	for (TIndexIndividual i = 0; i < a_populationSize;i++){
		a_p_reproductionProbs[i] = ReproductionProba(i+1);
		}
	}

TIndexIndividual* CPopulation::SelectedFuturProgeny(){
	a_p_prng->multinomial((unsigned int*)perIndividualOffspring, a_p_reproductionProbs , a_populationPresent.nbElts, a_populationPresent.nbElts);
	return perIndividualOffspring;
	}

void CPopulation::SortByFitness(){
	qsort(&a_populationPresent.p_elt[0], a_populationPresent.nbElts, sizeof(CIndividual*) , CompareIndividualsByFitness);

	}

void CPopulation::checkAvailableOccupiedPositions(TIndexIndividual* perIndividualOffspring){
	TIndexIndividual indexAvailables = 0;
	TIndexIndividual indexOccupieds = 0;
	for (TIndexIndividual i = 0; i < a_populationSize;i++){
		if (perIndividualOffspring[i] == 0){
			a_availablePositions[indexAvailables]=i;
			indexAvailables++;
			}
		else{
			a_occupiedPositions[indexOccupieds]=i;
			indexOccupieds++;
			}
		}
	a_nbAvailablePositions= indexAvailables;
	a_nbOccupiedPositions= indexOccupieds;
	}

void CPopulation::GiveBirth (TIndexIndividual* perIndividualOffspring){
	checkAvailableOccupiedPositions(perIndividualOffspring);
	TIndexIndividual indexAvailables = 0;
	for (TIndexIndividual i = 0; i < a_nbOccupiedPositions;i++){
		for (TIndexIndividual j = 0; j < perIndividualOffspring[a_occupiedPositions[i]]-1;j++){
			a_populationPresent.p_elt[a_availablePositions[indexAvailables]]->CopyOtherIndividual(a_populationPresent.p_elt[a_occupiedPositions[i]]);
			indexAvailables++;
			}
		}
	}

void CPopulation::MutateNextGeneration(TSimulationParameters * simulationParameters){
	for (TIndexIndividual i = 0; i < a_populationSize;i++){
		a_populationPresent.p_elt[i]->DoAllLargeMutations(simulationParameters);
		a_populationPresent.p_elt[i]->DoAllPointMutations(simulationParameters);
		}
	}

void CPopulation::SortNewGeneration(){
	qsort(&a_populationPresent.p_elt[0], a_populationPresent.nbElts, sizeof(CIndividual*) , CompareIndividualsByPosition);
	}

void CPopulation::ShuffleGenomes(){
	for (TIndexIndividual i = 0; i < a_populationSize;i++){
		a_populationPresent.p_elt[i]->Shuffle_genome();
		}
	}

void CPopulation::ComputeNextGeneration(TSimulationParameters * simulationParameters){
	SortByFitness();
	//printf("Sorted by fitness\n");
	GiveBirth(SelectedFuturProgeny());
	//printf("given birth\n");
	MutateNextGeneration(simulationParameters);
	//printf("next generation mutated\n");
	if (simulationParameters->sortByPosition) SortNewGeneration();
	//printf("pop sorted \n");
	if (simulationParameters->shuffleGenome)  ShuffleGenomes();
	}

void CPopulation::ComputeIndividualFitness(TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData, CFitnessEvaluator * fitnessComputer,TBoolean saveBestIndividual){
	for (TIndexIndividual i = 0; i < a_populationSize;i++){
		a_populationPresent.p_elt[i]->ComputeFitness(dataStream, startObservationPosInData, stopObservationPosInData, fitnessComputer);
		}
	if(saveBestIndividual){
		a_bestOldIdividual->ComputeFitness(dataStream, startObservationPosInData, stopObservationPosInData, fitnessComputer);
		}
	}

TIndexIndividual CPopulation::GetBestIndividualIndex(){
  TIndexIndividual bestIndividual = -1;
  TFitness         bestFitness    = -FLT_MAX;
  for (TIndexIndividual i = 0; i < a_populationSize;i++){
    if (a_populationPresent.p_elt[i]->a_fitness >= bestFitness){
      bestIndividual = i;
      bestFitness    = a_populationPresent.p_elt[i]->a_fitness;
      }
	}
  return(bestIndividual);
}

void CPopulation::SaveOrUpdateBestIndividual(TBoolean saveBestIndividual){
  if(saveBestIndividual){
    CIndividual* bestNewIndividual = a_populationPresent.p_elt[GetBestIndividualIndex()];
    if (a_bestOldIdividual->a_fitness <= bestNewIndividual->a_fitness){
      a_bestOldIdividual->CopyOtherIndividual(bestNewIndividual);
      }
    else{
      a_populationPresent.p_elt[a_p_prng->uniform(0,a_populationSize)]->CopyOtherIndividual(a_bestOldIdividual);
    }
  }
}

void CPopulation::EvolveDuringTsteps(TSimulationParameters * simulationParameters,TGenerationsIndex T, TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData , CFitnessEvaluator * fitnessComputer ){
	for (TGenerationsIndex i = 0; i< T; i++){
		//printf("________________________ \n");
		ComputeNextGeneration(simulationParameters);
		//printf("next generation computed \n");
		ComputeIndividualFitness(dataStream, startObservationPosInData, stopObservationPosInData , fitnessComputer, simulationParameters->saveBestIndividual);
    	//printf("fitnesses computed \n");
    	SaveOrUpdateBestIndividual(simulationParameters->saveBestIndividual);
		//printf("best individual updated \n");
		a_paramCurentGenerationindex ++;
		//printf("next generation increased \n");
		}
	}

int CPopulation::CompareIndividualsByFitness(const void* A_i, const void* B_i){
	CIndividual ** A1 = (CIndividual **) A_i;
	CIndividual ** B1 = (CIndividual **) B_i;
	CIndividual * A = *A1 ;
	CIndividual * B = *B1 ;
	return (int)(A->a_fitness - B->a_fitness);
	}

int CPopulation::CompareIndividualsByPosition(const void* A_i, const void* B_i){
	CIndividual ** A1 = (CIndividual **) A_i;
	CIndividual ** B1 = (CIndividual **) B_i;
	CIndividual * A = *A1 ;
	CIndividual * B = *B1 ;
	return (int)(A->a_position - B->a_position);
	}


void CPopulation::SavePopulationFeatures(TFileName nameSaveFile){
	for (TIndexIndividual i = 0; i < a_populationSize;i++){
		a_populationPresent.p_elt[i]->SaveFeatures(nameSaveFile,a_paramCurentGenerationindex);
		}
	}

void CPopulation::SaveIndividualGenome(TFileName nameSaveFile,TIndexIndividual individual){
	a_populationPresent.p_elt[individual]->SaveGenome(nameSaveFile);
	}

void CPopulation::SaveIndividualFeatures(TFileName nameSaveFile,TIndexIndividual individual){
	a_populationPresent.p_elt[individual]->SaveFeatures(nameSaveFile, a_paramCurentGenerationindex);
	}

void CPopulation::SaveIndividualPhenotype( CFitnessEvaluator * fitnessComputer, TFileName nameSaveFile,TIndexIndividual individual){
	a_populationPresent.p_elt[individual]->SavePhenotype( fitnessComputer,  nameSaveFile);
	}

void CPopulation::SaveIndividualClassifiedData(TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData ,CFitnessEvaluator * fitnessComputer, TFileName nameSaveFile,TIndexIndividual individual){
	a_populationPresent.p_elt[individual]->SaveClassifiedData(dataStream, startObservationPosInData, stopObservationPosInData ,fitnessComputer, nameSaveFile);
}

void CPopulation::ModifyPseudogenesAllPopulation(TMutationLaw** newPseudogenePDF, TBoundaryConditions** geneBoundaryConditions){
	for (TIndexIndividual i = 0; i < a_populationSize;i++){
		a_populationPresent.p_elt[i]->ModifyPseudogenes(newPseudogenePDF, geneBoundaryConditions);
	}
}

void CPopulation::ModifyPseudogenesOneIndividual(TMutationLaw** newPseudogenePDF, TBoundaryConditions** geneBoundaryConditions,TIndexIndividual individual){
	a_populationPresent.p_elt[individual]->ModifyPseudogenes(newPseudogenePDF, geneBoundaryConditions);
}

void CPopulation::OrderGenomesAllPopulation(){
	for (TIndexIndividual i = 0; i < a_populationSize;i++){
		a_populationPresent.p_elt[i]->OrderGenome();
	}
}
void CPopulation::OrderGenomeOneIndividual(TIndexIndividual individual){
	a_populationPresent.p_elt[individual]->OrderGenome();
	}

void CPopulation::ShuffleGenomeAllPopulation(){
	for (TIndexIndividual i = 0; i < a_populationSize;i++){
		a_populationPresent.p_elt[i]->Shuffle_genome();
		}
	}
void CPopulation::ShuffleGenomeOneIndividual(TIndexIndividual individual){
	a_populationPresent.p_elt[individual]->Shuffle_genome();
	}

