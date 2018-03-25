#include <zlib.h>
#include <errno.h>
#include <string.h>
#include "CSimulation.h"
#include "CPopulation.h"
#include "CIndividual.h"




CSimulation::CSimulation(TSimulationParameters* parameters){
	a_parameters = parameters;
	a_simulationState = (TSimulationState*) malloc(sizeof(TSimulationState));
	a_p_prng = new CPrng(a_parameters->prngSeed);
	a_simulationState->r_population = NULL;
	a_startObservationPosInData = 0;
	a_stopObservationPosInData = 0;
	a_fitnessEvaluator = NULL;
	InitRandomPopulation();
	AllocateDataStreamMem();
	InitializeFitnessFunction();
	printf("parameters s=%f seed=%f u_Dup=%f u_Del=%f u_Inv=%f u_tra=%f u_subs=%f u_pIns=%f u_pDel=%f e=%i N=%i Gi=%i multiMut=%i S=%i\n",a_parameters->selectionPressure ,a_parameters->prngSeed ,a_parameters->largeDuplicationRate,a_parameters->largeDeletionRate,a_parameters->largeInvertionRate,a_parameters->largeTranslocationRate,a_parameters->pointSubstitutionRate, a_parameters->pointInsertionRate, a_parameters->pointDeletionRate ,a_parameters->saveBestIndividual,a_parameters->populationSize,a_parameters->numberOfInitialGenes,a_parameters->moreThanOnePointSubstitution,a_parameters->sizeDataBuffer);
}


CSimulation::CSimulation(gzFile* backup_file,FILE* backup_state_rng, TSimulationParameters* parameters){//aumentar parametro CPRNG
	a_parameters = parameters;
	a_simulationState = (TSimulationState*) malloc(sizeof(TSimulationState));
	gzread((gzFile) backup_file, &a_startObservationPosInData, sizeof(TIndexPoint));
	gzread((gzFile) backup_file, &a_stopObservationPosInData , sizeof(TIndexPoint));
	a_p_prng = new CPrng(a_parameters->prngSeed);
	a_p_prng->loadState(backup_state_rng);
	AllocateDataStreamMem();
	InitializeFitnessFunction();
	a_simulationState->r_population = new CPopulation(backup_file , a_p_prng);
	}

void CSimulation::save( gzFile* backup_file , FILE* backup_state_rng){
	gzwrite((gzFile) backup_file, &a_startObservationPosInData, sizeof(TIndexPoint));
	gzwrite((gzFile) backup_file, &a_stopObservationPosInData , sizeof(TIndexPoint));
	a_p_prng->saveState(backup_state_rng);
	a_simulationState->r_population->save(backup_file);
	}

void CSimulation::InitializeFitnessFunction(){
	switch(a_parameters->fitnessMode){
		case FITNESS_MODE_TEST_1:{
			a_fitnessEvaluator = new CFitnessModeTest1( a_parameters->fitnessMode,1, a_parameters->boundaryConditions[INDEX_DIM]->max  - a_parameters->boundaryConditions[INDEX_DIM]->min  );
			break;
		}
		case FITNESS_MODE_TEST_2:{
			a_fitnessEvaluator = new CFitnessModeTest1( a_parameters->fitnessMode,1, a_parameters->boundaryConditions[INDEX_DIM]->max  - a_parameters->boundaryConditions[INDEX_DIM]->min  );
			break;
		}
		case FITNESS_MODE_TEST_3:{
			a_fitnessEvaluator = new CFitnessModeTest1( a_parameters->fitnessMode,1, a_parameters->boundaryConditions[INDEX_DIM]->max  - a_parameters->boundaryConditions[INDEX_DIM]->min  );
			break;
		}
		case FITNESS_MODE_CLUSTERING:{
			a_fitnessEvaluator = new CFitnessCluster(a_p_prng,a_parameters->genomeSizeLimit,a_parameters->nonCodingGenomeFitness,a_parameters->normExponent,a_parameters->unknownDimentionRandomValuesGenerator, a_parameters->boundaryConditions[INDEX_CLUST]->max  - a_parameters->boundaryConditions[INDEX_CLUST]->min , a_parameters->boundaryConditions[INDEX_DIM]->max - a_parameters->boundaryConditions[INDEX_DIM]->min);
			break;
		}
		case FITNESS_MODE_CLUSTERING_KMEANS:{
			a_fitnessEvaluator = new CFitnessClusterKmeans(a_p_prng,a_parameters->genomeSizeLimit,a_parameters->nonCodingGenomeFitness,a_parameters->normExponent,a_parameters->unknownDimentionRandomValuesGenerator, a_parameters->kmeansIterations, a_parameters->boundaryConditions[INDEX_CLUST]->max  - a_parameters->boundaryConditions[INDEX_CLUST]->min  , a_parameters->boundaryConditions[INDEX_DIM]->max - a_parameters->boundaryConditions[INDEX_DIM]->min  , a_parameters->sizeDataBuffer);
			break;
		}
	}
}
void CSimulation::ReallocateDataStreamMem(){
	freeDataStreamMem(a_dataStream);
	AllocateDataStreamMem();
}
void CSimulation::AllocateDataStreamMem(){
	a_dataStream = (TData*)malloc(sizeof(TData));
	a_dataStream->p_arrayPoints = (TArrayPoints*)malloc(sizeof(TArrayPoints));
	//int nbOfPoints =1;
	//int nbCoord = 10;
	a_dataStream->p_arrayPoints->p_elt = (TPoint *)malloc(a_parameters->sizeDataBuffer*sizeof(TPoint));
	a_dataStream->p_arrayPoints->nbElts = a_parameters->sizeDataBuffer;

	for (int i = 0;i<a_dataStream->p_arrayPoints->nbElts;i++){
		a_dataStream->p_arrayPoints->p_elt[i].p_coordinates = (TPointCoordinatesArray *)malloc(sizeof(TPointCoordinatesArray));
		a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt = (TCoordinate * )malloc(a_parameters->sizeDataPoint*sizeof(TCoordinate));
		a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts=0;//must be changed dynamically when non void data is added to the data stream
		a_dataStream->p_arrayPoints->p_elt[i].label   = -1;
		a_dataStream->p_arrayPoints->p_elt[i].cluster = -1;
		for (int j = 0; j<a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++){
			a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE] = 0;
			a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE] = 0;
			}
		}
	}

void CSimulation::LoadParameters(TSimulationParameters* parameters){
	a_parameters = parameters ;
	if (a_simulationState->r_population != NULL){
		ChangePopulationSize(parameters->populationSize );
		ChangeSelectionPressure( parameters->selectionPressure);
		}
	}

void CSimulation::freeMutationLawsArray(TMutationLaw** mutationLaws){
	if (mutationLaws!=NULL){
		for (TGeneElement i = 0; i < a_parameters->geneSize; i++){
			if (mutationLaws[i] != NULL) freeMutationLaw(mutationLaws[i]);
			}
		free(mutationLaws);
		}
	}

void CSimulation::freeIntergenicCutParameter(TIntergenicCut* intergenicCutParameter){
	if (intergenicCutParameter!=NULL){
		free(intergenicCutParameter->genesAbstractOrder);
		free(intergenicCutParameter);
	}
}

void CSimulation::freeMutationLaw(TMutationLaw* mutationLaw){
	if (mutationLaw->mutationLaw->p_transitionMatrix->matrixSize > 0){
		for (TIndexTransitionProbability i = 0;i < mutationLaw->mutationLaw->p_transitionMatrix->matrixSize;i++){
			free(mutationLaw->mutationLaw->p_transitionMatrix->matrix[i]);
			}
		free(mutationLaw->mutationLaw->p_transitionMatrix->matrix);
		}
	free(mutationLaw->mutationLaw->p_transitionMatrix);
	free(mutationLaw->mutationLaw);
	free(mutationLaw);
	}

void CSimulation::freeBoundaryConditions(TBoundaryConditions** boundaryConditions){
	for (TIndexGeneElement i = 0; i<a_parameters->geneSize;i++){
		free(boundaryConditions[i]);
		}
	free(boundaryConditions);
	}

void CSimulation::freeGeneEleMutationProbLaw(TGeneEleMutationProbsLaw* geneEleMutationProbsLaw){
	free(geneEleMutationProbsLaw->probGeneElementMutation);
	free(geneEleMutationProbsLaw->probGeneElementMutationWheelOfFortune);
	free(geneEleMutationProbsLaw);
	}


void CSimulation::freeDataStreamMem(TData* dataStream){
	for (int i = 0;i<a_dataStream->p_arrayPoints->nbElts;i++){
		free(dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt);
		free(a_dataStream->p_arrayPoints->p_elt[i].p_coordinates);
		}
	free(dataStream->p_arrayPoints->p_elt);
	free(dataStream->p_arrayPoints);
	free(dataStream);
	}



CSimulation::~CSimulation(){

	freeDataStreamMem(a_dataStream);
	freeMutationLawsArray(a_parameters->mutationLaws);
	freeMutationLawsArray(a_parameters->initPositionsPDF);
	freeBoundaryConditions(a_parameters->boundaryConditions);
	freeGeneEleMutationProbLaw(a_parameters->geneEleMutationProbsLaw);
	freeIntergenicCutParameter(a_parameters->intergenicCutParameter);
	if (a_simulationState->r_population!=NULL) free(a_simulationState->r_population);
	free(a_simulationState);
	free(a_parameters);
	freeMutationLaw(a_parameters->unknownDimentionRandomValuesGenerator);
	if (a_fitnessEvaluator != NULL) delete a_fitnessEvaluator;
	delete a_p_prng;
}

void CSimulation::InitRandomPopulation(){
	if (a_simulationState->r_population == NULL)a_simulationState->r_population = new CPopulation(a_parameters->populationSize, a_parameters->selectionPressure, a_p_prng, a_parameters->currentGenerationIndex);
	a_simulationState->r_population->GenerateRandomPopulation(a_parameters->numberOfInitialGenes, a_parameters->genomeSizeLimit,a_parameters->geneSize,a_parameters->initPositionsPDF, a_parameters->boundaryConditions);
	}

void CSimulation::PrintGenomes(){
	a_simulationState->r_population->PrintAllGenomes();
	}

void CSimulation::ChangePopulationSize(TPopulationSize newPopulationSize){
	TPopulationSize oldPopulationSize = a_simulationState->r_population->a_populationSize;
	a_simulationState->r_population->a_populationSize = newPopulationSize;
	a_simulationState->r_population->a_populationPresent.p_elt = (CIndividual ** ) realloc(a_simulationState->r_population->a_populationPresent.p_elt, newPopulationSize*sizeof(CIndividual*));
	a_simulationState->r_population->a_populationPresent.nbElts = newPopulationSize;
	a_simulationState->r_population->perIndividualOffspring = (TIndexIndividual*) realloc(a_simulationState->r_population->perIndividualOffspring, newPopulationSize * sizeof(TIndexIndividual));
	a_simulationState->r_population->a_availablePositions = (TIndexIndividual*) realloc(a_simulationState->r_population->a_availablePositions, newPopulationSize * sizeof(TIndexIndividual));
	a_simulationState->r_population->a_occupiedPositions = (TIndexIndividual*) realloc(a_simulationState->r_population->a_occupiedPositions, newPopulationSize * sizeof(TIndexIndividual));
	a_simulationState->r_population->a_p_reproductionProbs=(TReproductionProb*)realloc(a_simulationState->r_population->a_p_reproductionProbs, newPopulationSize*sizeof(TReproductionProb));
	a_simulationState->r_population->GenerateReproductionProbas();
	a_simulationState->r_population->GenerateRandomPopulation(a_parameters->numberOfInitialGenes,a_parameters->genomeSizeLimit,a_parameters->geneSize,a_parameters->initPositionsPDF, a_parameters->boundaryConditions,oldPopulationSize, newPopulationSize);
}

void CSimulation::ChangeSelectionPressure( TSelectionStrength newSelectionPressure){
	a_simulationState->r_population->a_selectionPressure =newSelectionPressure;
	a_simulationState->r_population->GenerateReproductionProbas();
}

TGenome* CSimulation::GetGenomeBestIndivLastGeneration(){
	a_simulationState->r_population->SortByFitness();
	return a_simulationState->r_population->a_populationPresent.p_elt[a_simulationState->r_population->a_populationPresent.nbElts-1]->a_genome;
	}

void CSimulation::EvolveDuringTsteps(TGenerationsIndex T){
	a_simulationState->r_population->EvolveDuringTsteps(a_parameters, T, a_dataStream, a_startObservationPosInData, a_stopObservationPosInData, (CFitnessEvaluator*) a_fitnessEvaluator);
	a_parameters->currentGenerationIndex+=T;
	}

TIndexPoint CSimulation::SetPoint(TIndexCoordinate pointDimentionality){
	//if the database is empty
	if (a_stopObservationPosInData == 0){
		a_dataStream->p_arrayPoints->p_elt[a_stopObservationPosInData].p_coordinates->nbElts = pointDimentionality;
		a_stopObservationPosInData++;
		return 0;
		}
	//if all the memory is full and start<stop
	if (a_stopObservationPosInData == a_parameters->sizeDataBuffer){
		a_dataStream->p_arrayPoints->p_elt[0].p_coordinates->nbElts = pointDimentionality;
		a_stopObservationPosInData = 1;
		a_startObservationPosInData = 1;
		return 0;
		}
	//if we didnt fill the entire dataset memory space
	if (a_startObservationPosInData < a_stopObservationPosInData){
		a_dataStream->p_arrayPoints->p_elt[a_stopObservationPosInData].p_coordinates->nbElts = pointDimentionality;
		a_stopObservationPosInData++;
		return a_stopObservationPosInData-1;
		}
	//if all the memory is full and start=stop=end
	if (a_startObservationPosInData == a_parameters->sizeDataBuffer-1){
		a_startObservationPosInData = 0;
		a_stopObservationPosInData ++;
		return a_stopObservationPosInData-1;
		}
	//if all the memory is full and start = stop
	if (a_stopObservationPosInData == a_startObservationPosInData){
		a_dataStream->p_arrayPoints->p_elt[a_stopObservationPosInData].p_coordinates->nbElts = pointDimentionality;
		a_stopObservationPosInData++;
		a_startObservationPosInData++;
		return a_stopObservationPosInData-1;
		}
	printf("errror-------------------------------------\n");
	return(-1);
	}

void CSimulation::SetPointCoordinate(TCoordinateElement* coordinates, TIndexCoordinate coordinateIndexInPoint, TIndexPoint pointIndex ){

	for (TIndexCoordinateElement j = 0; j<COORDINATE_LENGTH; j++){
			a_dataStream->p_arrayPoints->p_elt[pointIndex].p_coordinates->p_elt[coordinateIndexInPoint].coordinate[j] = coordinates[j];
		}
	}

void  CSimulation::PrintDataSet(){
	if  (a_startObservationPosInData < a_stopObservationPosInData){
		for (TIndexPoint i = a_startObservationPosInData; i<a_stopObservationPosInData;i++){
			for(TIndexCoordinate j = 0; j<a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				printf("%d %d\n",a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
				}
			}
		}
	if (a_startObservationPosInData == a_stopObservationPosInData){

		for (TIndexPoint i = 0; i<a_stopObservationPosInData;i++){
			for(TIndexCoordinate j = 0; j<=a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				printf("%d %d\n",a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
				}
			}

		for (TIndexPoint i = a_startObservationPosInData; i<a_dataStream->p_arrayPoints->nbElts;i++){
			for(TIndexCoordinate j = 0; j<=a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++ ){
				printf("%d %d\n",a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE]);
				}
			}
		}
	}

void CSimulation::ChangeSeed(Tseed  newSeed){
	//printf("new seed %lg\n",newSeed);
	a_parameters->prngSeed = newSeed;
	a_p_prng = new CPrng(a_parameters->prngSeed);
	a_simulationState->r_population->SetPRNG(a_p_prng);
	if (a_fitnessEvaluator != NULL) a_fitnessEvaluator->a_prng = a_p_prng;
	}

void CSimulation::SaveIndividualGenome(TFileName nameSaveFile,TIndexIndividual individual){
	a_simulationState->r_population->SaveIndividualGenome( nameSaveFile, individual);
	}

void CSimulation::SaveIndividualFeatures(TFileName nameSaveFile, TIndexIndividual individual){
	a_simulationState->r_population->SaveIndividualFeatures(nameSaveFile,individual);
	}

void CSimulation::SaveIndividualPhenotype( TFileName nameSaveFile,TIndexIndividual individual){
	a_simulationState->r_population->SaveIndividualPhenotype(a_fitnessEvaluator,nameSaveFile,individual);
	}

void CSimulation::SavePopulationFeatures(TFileName nameSaveFile){
	a_simulationState->r_population->SavePopulationFeatures(nameSaveFile);
	}

void CSimulation::SaveIndividualClassifiedData(TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData , TFileName nameSaveFile,TIndexIndividual individual){
  a_simulationState->r_population->SaveIndividualClassifiedData(dataStream,  startObservationPosInData,  stopObservationPosInData ,  a_fitnessEvaluator, nameSaveFile, individual);
	}

void CSimulation::ModifyPseudogenesAllPopulation(TMutationLaw** newPseudogenePDF){
	a_simulationState->r_population->ModifyPseudogenesAllPopulation(newPseudogenePDF, a_parameters->boundaryConditions);
	freeMutationLawsArray(newPseudogenePDF);
	}

void CSimulation::ModifyPseudogenesOneIndividual(TMutationLaw** newPseudogenePDF, TIndexIndividual individual){
	a_simulationState->r_population->ModifyPseudogenesOneIndividual(newPseudogenePDF, a_parameters->boundaryConditions, individual);
	freeMutationLawsArray(newPseudogenePDF);
	}

void CSimulation::OrderGenomesAllPopulation(){
	a_simulationState->r_population->OrderGenomesAllPopulation();
	}

void CSimulation::OrderGenomeOneIndividual(TIndexIndividual individual){
	a_simulationState->r_population->OrderGenomeOneIndividual(individual);
	}

void CSimulation::ShuffleGenomeAllPopulation(){
	a_simulationState->r_population->ShuffleGenomeAllPopulation();
	}
void CSimulation::ShuffleGenomeOneIndividual(TIndexIndividual individual){
	a_simulationState->r_population->ShuffleGenomeOneIndividual( individual);
	}
