#include "CIndividual.h"
#include "CPopulation.h"
#include "CSimulation.h"
#include "evoevo_def_state_pop.h"
#include "evoevo_def_base.h"
#include "evoevo_def_parameters.h"
#include "CPrng.h"
#include "evoevo_def_PyC.h"
#ifdef __cplusplus
extern "C" {
#endif
#include <Python.h>//include the "Python.h" header before any other include
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <zlib.h>
#include <errno.h>


static CSimulation* SimulationPythonToC(PyObject* args){
	CSimulation* simulation;
	PyObject* capsule;
	if (!PyArg_ParseTuple(args, "O", &capsule)){
		return NULL;
	}
	simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	return simulation;
}

void SimulationCapsuleDestructor(PyObject* capsule){
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	delete simulation;
}

static void AllocateSimulationParameter(TSimulationParameters* par ){
	TMutationLaw ** mutlaws = (TMutationLaw **)  malloc( par->geneSize *sizeof(TMutationLaw*));
	TMutationLaw ** initialDistributionLawPoints = (TMutationLaw **)  malloc( par->geneSize *sizeof(TMutationLaw*));
	TBoundaryConditions**  boundaryConditions = (TBoundaryConditions**)  malloc( par->geneSize *sizeof(TBoundaryConditions*));
	
	TGeneEleMutationProbsLaw* geneEleMutationProbsLaw = (TGeneEleMutationProbsLaw*) malloc(sizeof(TGeneEleMutationProbsLaw));
	geneEleMutationProbsLaw->probGeneElementMutation = (TMutationRate*) malloc(par->geneSize*sizeof(TMutationRate));
	geneEleMutationProbsLaw->probGeneElementMutationWheelOfFortune = (TMutationRate*) malloc((par->geneSize+1)*sizeof(TMutationRate));
	
	for (TIndexGeneElement i = 0; i<par->geneSize;i++){
		mutlaws[i] = (TMutationLaw *)  malloc(sizeof(TMutationLaw));
		mutlaws[i]->mutationLaw = (TMutationLawParameters*)  malloc(sizeof(TMutationLawParameters));
		mutlaws[i]->mutationLaw->p_transitionMatrix = (TTransitionMatrix*) malloc(sizeof(TTransitionMatrix));
	}

	for (TIndexGeneElement i = 0; i<par->geneSize;i++){
		initialDistributionLawPoints[i] = (TMutationLaw *)  malloc(sizeof(TMutationLaw));
		initialDistributionLawPoints[i]->mutationLaw = (TMutationLawParameters*)  malloc(sizeof(TMutationLawParameters));
		initialDistributionLawPoints[i]->mutationLaw->p_transitionMatrix = (TTransitionMatrix*) malloc(sizeof(TTransitionMatrix));
	}

	for (TIndexGeneElement i = 0; i<par->geneSize;i++){
		boundaryConditions[i] = (TBoundaryConditions*)  malloc(sizeof(TBoundaryConditions));
	}
		
	par->intergenicCutParameter = (TIntergenicCut*) malloc(sizeof(TIntergenicCut));
	par->intergenicCutParameter->genesAbstractOrder = (TGene*) malloc(sizeof(TGene));
	
    par->unknownDimentionRandomValuesGenerator = (TMutationLaw *)  malloc(sizeof(TMutationLaw));
    par->unknownDimentionRandomValuesGenerator->mutationLaw = (TMutationLawParameters*)  malloc(sizeof(TMutationLawParameters));
    par->unknownDimentionRandomValuesGenerator->mutationLaw->p_transitionMatrix = (TTransitionMatrix*) malloc(sizeof(TTransitionMatrix));

	par	->boundaryConditions = boundaryConditions;
	par->mutationLaws = mutlaws;
	par->initPositionsPDF = initialDistributionLawPoints;
	par->geneEleMutationProbsLaw = geneEleMutationProbsLaw;
}

static void LoadSimpleParameters(gzFile* backup_file,TSimulationParameters* par){
	gzread((gzFile) backup_file, &par->geneSize, sizeof(TIndexGeneElement));
	gzread((gzFile) backup_file, &par->numberOfInitialGenes, sizeof(TIndexGene));
	gzread((gzFile) backup_file, &par->sizeDataBuffer, sizeof(TIndexPoint));
	gzread((gzFile) backup_file, &par->sizeDataPoint, sizeof(TIndexCoordinate));
	gzread((gzFile) backup_file, &par->largeDuplicationRate, sizeof(TMutationRate));
	gzread((gzFile) backup_file, &par->largeDeletionRate, sizeof(TMutationRate));
	gzread((gzFile) backup_file, &par->largeTranslocationRate, sizeof(TMutationRate));
    gzread((gzFile)backup_file, &par->largeInvertionRate, sizeof(TMutationRate));
    gzread((gzFile)backup_file, &par->pointSubstitutionRate, sizeof(TMutationRate));
	gzread((gzFile) backup_file, &par->pointInsertionRate, sizeof(TMutationRate));
    gzread((gzFile) backup_file, &par->pointDeletionRate, sizeof(TMutationRate));
    gzread((gzFile) backup_file, &par->moreThanOnePointSubstitution, sizeof( TBoolean));
    gzread((gzFile) backup_file, &par->populationSize, sizeof( TPopulationSize));
	gzread((gzFile) backup_file, &par->currentGenerationIndex, sizeof( TGenerationsIndex));
	gzread((gzFile) backup_file, &par->selectionPressure, sizeof( TSelectionStrength));
	gzread((gzFile) backup_file, &par->kmeansIterations, sizeof(TKmeansIterations));
	gzread((gzFile) backup_file, &par->normExponent, sizeof(TNormExponent));
	gzread((gzFile) backup_file, &par->fitnessMode, sizeof(TFitnessMode));
	gzread((gzFile) backup_file, &par->logFile, sizeof(TLogFile ));
	gzread((gzFile) backup_file, &par->prngSeed, sizeof(Tseed ));
	gzread((gzFile) backup_file, &par->genomeSizeLimit, sizeof(TIndexGene));
	gzread((gzFile) backup_file, &par->nonCodingGenomeFitness, sizeof(TFitness));
	gzread((gzFile) backup_file, &par->shuffleGenome,sizeof(TBoolean));
}

static void SaveSimpleParameters(gzFile* backup_file,TSimulationParameters* par){
	gzwrite((gzFile) backup_file, &par->geneSize, sizeof(TIndexGeneElement));
	gzwrite((gzFile) backup_file, &par->numberOfInitialGenes, sizeof(TIndexGene));
	gzwrite((gzFile) backup_file, &par->sizeDataBuffer, sizeof(TIndexPoint));
	gzwrite((gzFile) backup_file, &par->sizeDataPoint, sizeof(TIndexCoordinate));
	gzwrite((gzFile) backup_file, &par->largeDuplicationRate, sizeof(TMutationRate));
	gzwrite((gzFile) backup_file, &par->largeDeletionRate, sizeof(TMutationRate));
	gzwrite((gzFile) backup_file, &par->largeTranslocationRate, sizeof(TMutationRate));
    gzwrite((gzFile) backup_file, &par->largeInvertionRate, sizeof(TMutationRate));
    gzwrite((gzFile) backup_file, &par->pointSubstitutionRate, sizeof(TMutationRate));
	gzwrite((gzFile) backup_file, &par->pointInsertionRate, sizeof(TMutationRate));
    gzwrite((gzFile) backup_file, &par->pointDeletionRate, sizeof(TMutationRate));
    gzwrite((gzFile) backup_file, &par->moreThanOnePointSubstitution, sizeof( TBoolean));
    gzwrite((gzFile) backup_file, &par->populationSize, sizeof( TPopulationSize));
	gzwrite((gzFile) backup_file, &par->currentGenerationIndex, sizeof( TGenerationsIndex));
	gzwrite((gzFile) backup_file, &par->selectionPressure, sizeof( TSelectionStrength));
	gzwrite((gzFile) backup_file, &par->kmeansIterations, sizeof(TKmeansIterations));
	gzwrite((gzFile) backup_file, &par->normExponent, sizeof(TNormExponent));
	gzwrite((gzFile) backup_file, &par->sortByPosition, sizeof( TBoolean));
	gzwrite((gzFile) backup_file, &par->fitnessMode, sizeof(TFitnessMode));
	gzwrite((gzFile) backup_file, &par->logFile, sizeof(TLogFile ));
	gzwrite((gzFile) backup_file, &par->prngSeed, sizeof(Tseed ));
	gzwrite((gzFile) backup_file, &par->genomeSizeLimit, sizeof(TIndexGene));
	gzwrite((gzFile) backup_file, &par->nonCodingGenomeFitness, sizeof(TFitness));
	gzwrite((gzFile) backup_file, &par->shuffleGenome,sizeof(TBoolean));
}

static void LoadTransitionMatrix(gzFile* backup_file, TTransitionMatrix* tmatrix){
	gzread((gzFile) backup_file, &tmatrix->matrixSize, sizeof(TIndexTransitionProbability));
	tmatrix->matrix = (TTransitionProbability**) malloc(tmatrix->matrixSize*sizeof(TTransitionProbability*));
	for (TIndexTransitionProbability i = 0; i < tmatrix->matrixSize; i++) {
		tmatrix->matrix[i] = (TTransitionProbability*) malloc(tmatrix->matrixSize*sizeof(TTransitionProbability));
		for (TIndexTransitionProbability j = 0; j < tmatrix->matrixSize; j++){
			gzread((gzFile) backup_file, &tmatrix->matrix[i][j],sizeof(TTransitionProbability));
		}
	}
}

static void LoadMutationLawParameters(gzFile* backup_file, TMutationLawParameters* mutationLawParameters){
	gzread((gzFile) backup_file, &mutationLawParameters->min, sizeof(TGeneElement));
	gzread((gzFile) backup_file, &mutationLawParameters->max, sizeof(TGeneElement));
	gzread((gzFile) backup_file, &mutationLawParameters->mean, sizeof(TStatisticMomentum ));
	gzread((gzFile) backup_file, &mutationLawParameters->standardDeviation, sizeof(TStatisticMomentum ));
	gzread((gzFile) backup_file, &mutationLawParameters->law, sizeof(TLawFamilly));
	LoadTransitionMatrix(backup_file,mutationLawParameters->p_transitionMatrix);
}

static void SaveTransitionMatrix(gzFile* backup_file, TTransitionMatrix* tmatrix){
	gzwrite((gzFile) backup_file, &tmatrix->matrixSize, sizeof(TIndexTransitionProbability));
	for (TIndexTransitionProbability i = 0; i < tmatrix->matrixSize; i++) {
		for (TIndexTransitionProbability j = 0; j < tmatrix->matrixSize; j++){
			gzwrite((gzFile) backup_file, &tmatrix->matrix[i][j],sizeof(TTransitionProbability));
		}
	}
}

static void SaveMutationLawParameters(gzFile* backup_file, TMutationLawParameters* mutationLawParameters){
	gzwrite((gzFile) backup_file, &mutationLawParameters->min, sizeof(TGeneElement));
	gzwrite((gzFile) backup_file, &mutationLawParameters->max, sizeof(TGeneElement));
	gzwrite((gzFile) backup_file, &mutationLawParameters->mean, sizeof(TStatisticMomentum ));
	gzwrite((gzFile) backup_file, &mutationLawParameters->standardDeviation, sizeof(TStatisticMomentum ));
	gzwrite((gzFile) backup_file, &mutationLawParameters->law, sizeof(TLawFamilly));
	SaveTransitionMatrix(backup_file,mutationLawParameters->p_transitionMatrix);
}

static void SaveIntergenicCutParameter(gzFile* backup_file, TIntergenicCut* intergenicCutParameter,TIndexGeneElement size){
	gzwrite((gzFile) backup_file, &intergenicCutParameter->cut, sizeof(TBoolean));
	gzwrite((gzFile) backup_file, &intergenicCutParameter->mincutpoint, sizeof(TGeneElement));
	gzwrite((gzFile) backup_file, &intergenicCutParameter->maxcutpoint, sizeof(TGeneElement));
	for(TGeneElement i = 0; i<size; i++){
		gzwrite((gzFile) backup_file, &intergenicCutParameter->genesAbstractOrder[i], sizeof(TGeneElement));
	}
}

static void  LoadIntergenicCutParameter(gzFile* backup_file, TIntergenicCut* intergenicCutParameter,TIndexGeneElement size){
	gzread((gzFile) backup_file, &intergenicCutParameter->cut, sizeof(TBoolean));
	gzread((gzFile) backup_file, &intergenicCutParameter->mincutpoint, sizeof(TGeneElement));
	gzread((gzFile) backup_file, &intergenicCutParameter->maxcutpoint, sizeof(TGeneElement));
	for(TGeneElement i = 0; i<size; i++){
		gzread((gzFile) backup_file, &intergenicCutParameter->genesAbstractOrder[i], sizeof(TGeneElement));
	}
}

static void LoadAllMutationLawParameters(gzFile* backup_file, TMutationLaw** mutationLaws, TIndexGeneElement size){
	for (TIndexGeneElement i = 0;i < size; i++){
		LoadMutationLawParameters(backup_file, mutationLaws[i]->mutationLaw);
	}
}

static void SaveAllMutationLawParameters(gzFile* backup_file, TMutationLaw** mutationLaws, TIndexGeneElement size){
	for (TIndexGeneElement i = 0; i < size; i++){
		SaveMutationLawParameters(backup_file, mutationLaws[i]->mutationLaw);
		}
	}

static void LoadGeneBoundaryConditions(gzFile* backup_file, TBoundaryConditions* boundaryConditions){
	gzread((gzFile) backup_file, &boundaryConditions->min, sizeof(TGeneElement));
	gzread((gzFile) backup_file, &boundaryConditions->max, sizeof(TGeneElement));
	}

static void SaveGeneBoundaryConditions(gzFile* backup_file, TBoundaryConditions* boundaryConditions){
	gzwrite((gzFile) backup_file, &boundaryConditions->min, sizeof(TGeneElement));
	gzwrite((gzFile) backup_file, &boundaryConditions->max, sizeof(TGeneElement));
	}

static void LoadAllGeneBoundaryConditions(gzFile* backup_file, TBoundaryConditions** boundaryConditions, TIndexGeneElement size){
	for (TIndexGeneElement i = 0; i < size; i++){
		LoadGeneBoundaryConditions( backup_file, boundaryConditions[i]);
		}
	}

static void SaveAllGeneBoundaryConditions(gzFile* backup_file, TBoundaryConditions** boundaryConditions, TIndexGeneElement size){
	for (TIndexGeneElement i = 0; i < size; i++){
		SaveGeneBoundaryConditions(backup_file, boundaryConditions[i]);
		}
	}

static void LoadAllGeneEleMutationProbsLaw(gzFile* backup_file, TGeneEleMutationProbsLaw* geneEleMutationProbsLaw, TIndexGeneElement size){
	for (TIndexGeneElement i = 0; i < size; i++){
		gzread((gzFile) backup_file, &geneEleMutationProbsLaw->probGeneElementMutation[i], sizeof(TMutationRate));
		}
	for (TIndexGeneElement i = 0; i < size+1; i++){
		gzread((gzFile) backup_file, &geneEleMutationProbsLaw->probGeneElementMutationWheelOfFortune[i], sizeof(TMutationRate));
		}
	}

static void SaveAllGeneEleMutationProbsLaw(gzFile* backup_file, TGeneEleMutationProbsLaw* geneEleMutationProbsLaw, TIndexGeneElement size){
	for (TIndexGeneElement i = 0; i < size; i++){
		gzwrite((gzFile) backup_file, &geneEleMutationProbsLaw->probGeneElementMutation[i], sizeof(TMutationRate));
		}
	for (TIndexGeneElement i = 0; i < size+1; i++){
		gzwrite((gzFile) backup_file, &geneEleMutationProbsLaw->probGeneElementMutationWheelOfFortune[i], sizeof(TMutationRate));
		}
	}

static TSimulationParameters* LoadParameters(TBackupFileName backup_file_name){
	gzFile * backup_file = (gzFile *) gzopen (backup_file_name, "r");
	if (! backup_file) {
		fprintf (stderr, "gzopen of '%s' failed: %s.\n", backup_file_name,
        strerror (errno));
        exit (EXIT_FAILURE);
		}
	TSimulationParameters* par = (TSimulationParameters*) malloc(sizeof(TSimulationParameters));
	LoadSimpleParameters(backup_file,par);
	AllocateSimulationParameter(par);
	LoadAllMutationLawParameters(backup_file, par->mutationLaws, par->geneSize);
	LoadAllMutationLawParameters(backup_file, par->initPositionsPDF, par->geneSize);
	LoadAllGeneBoundaryConditions(backup_file, par->boundaryConditions, par->geneSize);
	LoadAllGeneEleMutationProbsLaw(backup_file, par->geneEleMutationProbsLaw, par->geneSize);
	LoadMutationLawParameters(backup_file, par->unknownDimentionRandomValuesGenerator->mutationLaw);
	LoadIntergenicCutParameter(backup_file, par->intergenicCutParameter,par->geneSize);
	gzclose ((gzFile)backup_file);
	return par;
	}


static PyObject* InitSimulationFromBackup(PyObject* self, PyObject* args){
	TBackupFileName backup_file_name;
	TBackupFileName backup_state_rng_name;
	TBackupFileName params_backup_file_name;
	PyObject* capsule;
	if (!PyArg_ParseTuple(args, "Osss",&capsule,&backup_file_name,&backup_state_rng_name, &params_backup_file_name)){
		return NULL;
		}
	gzFile * backup_file = (gzFile *)gzopen (backup_file_name, "r");
	FILE* backup_state_rng = fopen ( backup_state_rng_name, "r");
	CSimulation* simulation = new CSimulation(backup_file,backup_state_rng,LoadParameters(params_backup_file_name));
	/*
	simulation->PrintGenomes();
	simulation->EvolveDuringTsteps(3);
	simulation->PrintGenomes();
	*/
	PyCapsule_SetPointer(capsule, simulation);
	/*
	printf("ssssssssss\n");
	CSimulation* simulation2 = (CSimulation*)PyCapsule_GetPointer(capsule,NULL);
	simulation2->PrintGenomes();
	*/
	Py_INCREF(Py_None);
	return Py_None;
	}


static void SaveParameters(TBackupFileName backup_file_name, TSimulationParameters* parameters){

		gzFile * backup_file =(gzFile *) gzopen (backup_file_name, "w");
		if (! backup_file) {
			fprintf (stderr, "gzopen of '%s' failed: %s.\n", backup_file_name,
            strerror (errno));
            exit (EXIT_FAILURE);
			}

		SaveSimpleParameters(backup_file,parameters);
		SaveAllMutationLawParameters(backup_file, parameters->mutationLaws, parameters->geneSize);
		SaveAllMutationLawParameters(backup_file, parameters->initPositionsPDF, parameters->geneSize);
		SaveAllGeneBoundaryConditions(backup_file, parameters->boundaryConditions, parameters->geneSize);
		SaveAllGeneEleMutationProbsLaw(backup_file, parameters->geneEleMutationProbsLaw, parameters->geneSize);
		SaveMutationLawParameters(backup_file, parameters->unknownDimentionRandomValuesGenerator->mutationLaw);
		SaveIntergenicCutParameter(backup_file, parameters->intergenicCutParameter,parameters->geneSize);
		gzclose ((gzFile)backup_file);
		}

static PyObject* SaveSimulation(PyObject* self, PyObject* args){
	TBackupFileName backup_file_name;
	TBackupFileName backup_state_rng_name;
	TBackupFileName params_backup_file_name;
	PyObject* capsule;
	if (!PyArg_ParseTuple(args, "sssO",&backup_file_name,&backup_state_rng_name, &params_backup_file_name,&capsule)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	SaveParameters(params_backup_file_name,simulation->a_parameters );
	gzFile * backup_file = (gzFile *) gzopen (backup_file_name, "w");
	FILE* backup_state_rng = fopen ( backup_state_rng_name, "w");
	simulation->save(backup_file , backup_state_rng);
	fclose(backup_state_rng);
	gzclose((gzFile)backup_file);
	Py_INCREF(Py_None);
	return Py_None;
	}


static void ParseBoundaryConditions(PyTupleObject* statisticalLawPython, TBoundaryConditions* boundaryConditions ){
	boundaryConditions->min = PyInt_AsLong(PyTuple_GetItem((PyObject*)statisticalLawPython, (Py_ssize_t)0));
	boundaryConditions->max = PyInt_AsLong(PyTuple_GetItem((PyObject*)statisticalLawPython, (Py_ssize_t)1));
	}
//*****************************************
static PyObject* BuildBoundaryConditions(TBoundaryConditions* boundaryConditions ){
	PyListObject* boundaryCondition = (PyListObject*)PyList_New(0);
	PyList_Append((PyObject*)boundaryCondition,PyInt_FromLong(boundaryConditions->min));
	PyList_Append((PyObject*)boundaryCondition,PyInt_FromLong(boundaryConditions->max));
	return (PyObject*)boundaryCondition;
	}

static void ParseMutationLawGeneElement (PyTupleObject* geneEleMutationProbsLawPython, TGeneEleMutationProbsLaw* geneEleMutationProbsLaw){

	PyListObject* geneElementMutationProbs = (PyListObject*)PyTuple_GetItem((PyObject*)geneEleMutationProbsLawPython, (Py_ssize_t)0);
	TIndexGeneElement sizeProbsArray = (TIndexGeneElement) PyList_Size((PyObject*) geneElementMutationProbs);
	PyListObject* geneElementWheelOfFortune = (PyListObject*)PyTuple_GetItem((PyObject*)geneEleMutationProbsLawPython, (Py_ssize_t)1);
	TIndexGeneElement sizeWheelOfFortune = (TIndexGeneElement) PyList_Size((PyObject*) geneElementWheelOfFortune);

	for (TIndexGeneElement i = 0; i < sizeProbsArray; i++ ){
		geneEleMutationProbsLaw->probGeneElementMutation[i] = PyFloat_AsDouble((PyObject*)PyList_GetItem((PyObject*)geneElementMutationProbs,(Py_ssize_t)i));
		}

	for (TIndexGeneElement i = 0; i < sizeWheelOfFortune; i++ ){
		geneEleMutationProbsLaw->probGeneElementMutationWheelOfFortune[i] = PyFloat_AsDouble((PyObject*)PyList_GetItem((PyObject*)geneElementWheelOfFortune,(Py_ssize_t)i));
		}
	}
//********************************************
static PyObject*  BuildMutationLawGeneElement (TGeneEleMutationProbsLaw* geneEleMutationProbsLaw,TIndexGeneElement size){
	PyListObject** geneEleMutationProbsLawPY = (PyListObject**)PyList_New(0);
	PyList_Append((PyObject*)geneEleMutationProbsLawPY, PyList_New(0));
	for (TIndexGeneElement i = 0; i < size; i++ ){
		PyList_Append((PyObject*)PyList_GetItem((PyObject*)geneEleMutationProbsLawPY,(Py_ssize_t)0), PyFloat_FromDouble(geneEleMutationProbsLaw->probGeneElementMutation[i]));
		}
	PyList_Append((PyObject*)geneEleMutationProbsLawPY, PyList_New(0));
	for (TIndexGeneElement i = 0; i < size+1; i++ ){
		PyList_Append((PyObject*)PyList_GetItem((PyObject*)geneEleMutationProbsLawPY,(Py_ssize_t)1), PyFloat_FromDouble(geneEleMutationProbsLaw->probGeneElementMutationWheelOfFortune[i]));
		}
	return (PyObject*)geneEleMutationProbsLawPY;
	}

static void ParseTransitionMatrix(PyListObject** transitionMatrixPython, TMutationLawParameters* statisticalLawC){
	TIndexTransitionProbability sizeMatrix = (TIndexTransitionProbability) PyList_Size((PyObject*)transitionMatrixPython);
	TIndexTransitionProbability size_loc;
	if (sizeMatrix > 0){
		statisticalLawC->p_transitionMatrix->matrixSize = sizeMatrix;
		statisticalLawC->p_transitionMatrix->matrix = (TTransitionProbability**) malloc(sizeMatrix*sizeof(TTransitionProbability*));
		for (TIndexTransitionProbability i = 0; i < sizeMatrix; i++) {
			size_loc = (TIndexTransitionProbability) PyList_Size((PyObject*)PyList_GetItem((PyObject*)transitionMatrixPython, (Py_ssize_t)i));
			statisticalLawC->p_transitionMatrix->matrix[i] = (TTransitionProbability*) malloc(size_loc*sizeof(TTransitionProbability));
			for (TIndexTransitionProbability j = 0; j < size_loc; j++){
				statisticalLawC->p_transitionMatrix->matrix[i][j] = PyFloat_AsDouble((PyObject*)PyList_GetItem((PyObject*)PyList_GetItem((PyObject*)transitionMatrixPython, (Py_ssize_t)i), (Py_ssize_t)j));
			}
		}
	}
	else {
		statisticalLawC->p_transitionMatrix->matrix = NULL;
		statisticalLawC->p_transitionMatrix->matrixSize = 0;
	}
}

//*****************************************
static PyObject* BuildTransitionMatrix(TMutationLawParameters* statisticalLawC){
	PyListObject** transitionMatrixPython = (PyListObject**)PyList_New(0);
	if (statisticalLawC->p_transitionMatrix->matrixSize > 0){
		for (TIndexTransitionProbability i = 0; i < statisticalLawC->p_transitionMatrix->matrixSize; i++) {
			PyList_Append((PyObject*)transitionMatrixPython, PyList_New(0));
			for (TIndexTransitionProbability j = 0; j < statisticalLawC->p_transitionMatrix->matrixSize; j++){
				PyList_Append((PyObject*)PyList_GetItem((PyObject*)transitionMatrixPython, (Py_ssize_t)i),PyFloat_FromDouble(statisticalLawC->p_transitionMatrix->matrix[i][j]));
				}
			}
		}
	return  (PyObject*)transitionMatrixPython;
	}

static void ParseStatisticalLaw(PyTupleObject* statisticalLawPython, TMutationLaw* statisticalLawC ){
	statisticalLawC->mutationLaw->min = PyInt_AsLong(PyTuple_GetItem((PyObject*)statisticalLawPython, (Py_ssize_t)0));
	statisticalLawC->mutationLaw->max = PyInt_AsLong(PyTuple_GetItem((PyObject*)statisticalLawPython, (Py_ssize_t)1));
	statisticalLawC->mutationLaw->mean = PyFloat_AsDouble(PyTuple_GetItem((PyObject*)statisticalLawPython, (Py_ssize_t)2));
	statisticalLawC->mutationLaw->standardDeviation = PyFloat_AsDouble(PyTuple_GetItem((PyObject*)statisticalLawPython, (Py_ssize_t)3));
	statisticalLawC->mutationLaw->law = PyInt_AsLong(PyTuple_GetItem((PyObject*)statisticalLawPython, (Py_ssize_t)4));
	ParseTransitionMatrix((PyListObject**)PyTuple_GetItem((PyObject*)statisticalLawPython, (Py_ssize_t)5), statisticalLawC->mutationLaw);
	}
//************************************
static  PyObject* BuildStatisticalLaw(TMutationLaw* statisticalLawC){
	PyListObject* mutLawPython = (PyListObject*)PyList_New(0);
	PyList_Append((PyObject*)mutLawPython,PyInt_FromLong(statisticalLawC->mutationLaw->min));//PyFloat_FromDouble(
	PyList_Append((PyObject*)mutLawPython,PyInt_FromLong(statisticalLawC->mutationLaw->max));
	PyList_Append((PyObject*)mutLawPython,PyInt_FromLong(statisticalLawC->mutationLaw->mean));
	PyList_Append((PyObject*)mutLawPython,PyInt_FromLong(statisticalLawC->mutationLaw->standardDeviation));
	PyList_Append((PyObject*)mutLawPython,PyInt_FromLong(statisticalLawC->mutationLaw->law));
	PyList_Append((PyObject*)mutLawPython,BuildTransitionMatrix(statisticalLawC->mutationLaw));
	return (PyObject*)mutLawPython;
	}

static void ParseAllGeneBoundaryConditions(PyTupleObject** statisticalLawsPython, TBoundaryConditions** boundaryConditionsC ){
	TIndexGeneElement size = (TIndexGeneElement) PyTuple_Size((PyObject*) statisticalLawsPython);
	for(TIndexGeneElement i = 0; i < size; i++){
		ParseBoundaryConditions((PyTupleObject*)PyTuple_GetItem((PyObject*)statisticalLawsPython, (Py_ssize_t)i),boundaryConditionsC[i]);
		}
	}
//****************************************
static PyObject* BuildAllGeneBoundaryConditions( TBoundaryConditions** boundaryConditionsC , TIndexGeneElement size){
	PyListObject** boundaryConditions = (PyListObject**)PyList_New(0);
	for(TIndexGeneElement i = 0; i < size; i++){
		PyList_Append((PyObject*)boundaryConditions,BuildBoundaryConditions(boundaryConditionsC[i]));
		}
	return (PyObject*)boundaryConditions;
	}

static void ParseAllGeneStatisticalLaws(PyTupleObject** statisticalLawsPython, TMutationLaw** statisticalLawsC ){
	TIndexGeneElement size = (TIndexGeneElement)PyTuple_Size((PyObject*) statisticalLawsPython);
	for(TIndexGeneElement i = 0; i < size; i++){
		ParseStatisticalLaw((PyTupleObject*)PyTuple_GetItem((PyObject*)statisticalLawsPython, (Py_ssize_t)i),statisticalLawsC[i]);
		}
	}
//*****************************************
static PyObject* BuildIntergenicCutParam( TIntergenicCut* intergenicCutParameter , TIndexGeneElement size){
	PyListObject* intergenicCutParamPython = (PyListObject*)PyList_New(0);
	PyList_Append((PyObject*)intergenicCutParamPython, PyInt_FromLong(intergenicCutParameter->cut));
	PyList_Append((PyObject*)intergenicCutParamPython, PyInt_FromLong(intergenicCutParameter->mincutpoint));
	PyList_Append((PyObject*)intergenicCutParamPython, PyInt_FromLong(intergenicCutParameter->maxcutpoint));
	for(TIndexGeneElement i = 0; i < size; i++){
		PyList_Append((PyObject*) intergenicCutParamPython,PyInt_FromLong(intergenicCutParameter->genesAbstractOrder->gene[i]));
		}
	return (PyObject*)intergenicCutParamPython;
	}

static void ParseIntergenicCutParam(PyTupleObject* intergenicCutParamPython, TIntergenicCut* intergenicCutParameter){
	TIndexGeneElement size = (TIndexGeneElement) PyTuple_Size((PyObject*)intergenicCutParamPython);
	//printf("%i\n",size);
	intergenicCutParameter->cut = (TBoolean)PyInt_AsLong(PyTuple_GetItem((PyObject*)intergenicCutParamPython, (Py_ssize_t)0));
	intergenicCutParameter->mincutpoint = PyInt_AsLong(PyTuple_GetItem((PyObject*)intergenicCutParamPython, (Py_ssize_t)1));
	intergenicCutParameter->maxcutpoint = PyInt_AsLong(PyTuple_GetItem((PyObject*)intergenicCutParamPython, (Py_ssize_t)2));
	//printf("%i %i %i\n",intergenicCutParameter->cut,intergenicCutParameter->mincutpoint,intergenicCutParameter->maxcutpoint);
	for (TIndexGeneElement i = 3; i < size; i++ ){
		intergenicCutParameter->genesAbstractOrder->gene[i-3] = PyInt_AsLong(PyTuple_GetItem((PyObject*)intergenicCutParamPython, (Py_ssize_t)i));
		//printf("%i\n",intergenicCutParameter->genesAbstractOrder->gene[i-3]);
		}
	}
//***************************************
static PyObject*  BuildAllGeneStatisticalLaws( TMutationLaw** statisticalLawsC , TIndexGeneElement size){
	PyListObject** mutLawsPython = (PyListObject**)PyList_New(0);
	for(TIndexGeneElement i = 0; i < size; i++){
		PyList_Append((PyObject*)mutLawsPython,BuildStatisticalLaw(statisticalLawsC[i]));
		}
	return (PyObject*)mutLawsPython;
	}

static PyObject*  Getparameters(PyObject* self,PyObject* args){
	CSimulation* simulation    = SimulationPythonToC(args);
	TSimulationParameters* par = simulation->a_parameters;
	BuildAllGeneStatisticalLaws(par->mutationLaws,par->geneSize);
	BuildAllGeneStatisticalLaws(par->initPositionsPDF,par->geneSize);
	BuildAllGeneBoundaryConditions(par->boundaryConditions,par->geneSize);
	BuildMutationLawGeneElement(par->geneEleMutationProbsLaw,par->geneSize);
	BuildStatisticalLaw(par->unknownDimentionRandomValuesGenerator);
	/*
	printf("blablablalbal\n");
	printf("%lg   dddddddddddd\n",par->prngSeed);
	printf("%i %i %i %i %i %i %f %f %f %f %f %f %f %i %i %i %f %i %f %i %s %f %i %i %f %f %f %i \n",par->geneSize, par->numberOfInitialGenes, par->sizeDataBuffer,par->sizeDataPoint,par->sizePopulationsMemory, par->genomeSizeLimit,par->largeDuplicationRate,par->largeDeletionRate,par->largeTranslocationRate,par->largeInvertionRate,par->pointSubstitutionRate,par->pointInsertionRate,par->pointDeletionRate,par->moreThanOnePointSubstitution,par->populationSize,par->currentGenerationIndex,par->selectionPressure,par->kmeansIterations,par->normExponent,par->fitnessMode,par->logFile,par->prngSeed,par->sortByPosition,par->permuteWhenRearrange,BuildAllGeneStatisticalLaws(par->mutationLaws,par->geneSize),par->nonCodingGenomeFitness ,par->dimentionNormWeight, par->residualDistanceWeight, par->shuffleGenome);

	printf("dsdfsfsdfsdfsdfsfs\n");
	*/
	return Py_BuildValue("iiiiifffffffiiififisfiOiOOOOfOi", par->geneSize,  par->numberOfInitialGenes,  par->sizeDataBuffer, par->sizeDataPoint, par->genomeSizeLimit, par->largeDuplicationRate, par->largeDeletionRate, par->largeTranslocationRate, par->largeInvertionRate, par->pointSubstitutionRate, par->pointInsertionRate, par->pointDeletionRate, par->moreThanOnePointSubstitution, par->populationSize, par->currentGenerationIndex, par->selectionPressure, par->kmeansIterations, par->normExponent, par->fitnessMode, par->logFile, par->prngSeed, par->sortByPosition, BuildIntergenicCutParam(par->intergenicCutParameter , par->geneSize), par->saveBestIndividual , BuildAllGeneStatisticalLaws(par->mutationLaws,par->geneSize), BuildAllGeneStatisticalLaws(par->initPositionsPDF,par->geneSize), BuildAllGeneBoundaryConditions(par->boundaryConditions,par->geneSize), BuildMutationLawGeneElement(par->geneEleMutationProbsLaw,par->geneSize),  par->nonCodingGenomeFitness , BuildStatisticalLaw(par->unknownDimentionRandomValuesGenerator),  par->shuffleGenome);
	//return Py_BuildValue("iiiiiifffffffiiififisfiiO",&par->geneSize, &par->numberOfInitialGenes, &par->sizeDataBuffer,&par->sizeDataPoint,&par->sizePopulationsMemory, &par->genomeSizeLimit,&par->largeDuplicationRate,&par->largeDeletionRate,&par->largeTranslocationRate,&par->largeInvertionRate,&par->pointSubstitutionRate,&par->pointInsertionRate,&par->pointDeletionRate,&par->moreThanOnePointSubstitution,&par->populationSize,&par->currentGenerationIndex,&par->selectionPressure,&par->kmeansIterations,&par->normExponent,&par->fitnessMode,&par->logFile,&par->prngSeed,&par->sortByPosition,&par->permuteWhenRearrange,BuildAllGeneStatisticalLaws(par->mutationLaws,par->geneSize));
	//return Py_BuildValue("iN",par->geneSize,BuildAllGeneStatisticalLaws(par->mutationLaws,par->geneSize));
	}

static void ParseSimulationParameters(PyObject* args, TSimulationParameters* par){

	PyTupleObject** mutLawsPython;
	PyTupleObject** initLawsPython;
	PyTupleObject** boundaryConditionsPython;
	PyTupleObject*  mutLawsGeneElementPython;
	PyTupleObject*  unknownDimentionRandomValuesGenerator;
	PyTupleObject*  intergenicCutParamPython;

	if (PyArg_ParseTuple(args, "iiiiifffffffiiififisfiOiOOOOfOi",&par->geneSize, &par->numberOfInitialGenes, &par->sizeDataBuffer,&par->sizeDataPoint, &par->genomeSizeLimit,&par->largeDuplicationRate,&par->largeDeletionRate,&par->largeTranslocationRate,&par->largeInvertionRate,&par->pointSubstitutionRate,&par->pointInsertionRate,&par->pointDeletionRate,&par->moreThanOnePointSubstitution,&par->populationSize,&par->currentGenerationIndex,&par->selectionPressure,&par->kmeansIterations,&par->normExponent,&par->fitnessMode,&par->logFile,&par->prngSeed,&par->sortByPosition,&intergenicCutParamPython,&par->saveBestIndividual ,&mutLawsPython,&initLawsPython,&boundaryConditionsPython,&mutLawsGeneElementPython, &par->nonCodingGenomeFitness ,&unknownDimentionRandomValuesGenerator, &par->shuffleGenome)){

		AllocateSimulationParameter(par);
		ParseAllGeneStatisticalLaws(mutLawsPython, par->mutationLaws);
		ParseAllGeneStatisticalLaws(initLawsPython, par->initPositionsPDF);
		ParseIntergenicCutParam(intergenicCutParamPython, par->intergenicCutParameter);
		ParseStatisticalLaw(unknownDimentionRandomValuesGenerator, par->unknownDimentionRandomValuesGenerator);

		ParseAllGeneBoundaryConditions(boundaryConditionsPython, par->boundaryConditions);
		ParseMutationLawGeneElement(mutLawsGeneElementPython, par->geneEleMutationProbsLaw);
		}
	}

static TSimulationParameters* InitSimulationParameter(PyObject* args){
	TSimulationParameters* par = (TSimulationParameters*) malloc(sizeof(TSimulationParameters));
	ParseSimulationParameters(args,par);
	return par;
	}

static PyObject* CreatSimulation_traductor(PyObject* self, PyObject* args){
	CSimulation* simulation = new CSimulation(InitSimulationParameter(args));
	PyObject* capsule = PyCapsule_New(simulation,NAME_CAPSULE, SimulationCapsuleDestructor);
	return capsule;
	}

static PyObject* SetParameters(PyObject* self, PyObject* args){
	PyObject* capsule;
	CSimulation* simulation;
	PyTupleObject** mutLawsPython;
	PyTupleObject** initLawsPython;
	PyTupleObject** boundaryConditionsPython;
	PyTupleObject*  mutLawsGeneElementPython;
	PyTupleObject*  unknownDimentionRandomValuesGenerator;
	TIndexGeneElement geneSize;
	TIndexGene numberOfInitialGenes;
	TIndexPoint sizeDataBuffer;
	TIndexCoordinate sizeDataPoint;
    TMutationRate largeDuplicationRate;
    TMutationRate largeDeletionRate;
    TMutationRate largeTranslocationRate;
    TMutationRate largeInvertionRate;
    TMutationRate pointSubstitutionRate;   
    TMutationRate pointInsertionRate;        
    TMutationRate pointDeletionRate;          
    TBoolean moreThanOnePointSubstitution; 
	TPopulationSize populationSize;
	TGenerationsIndex currentGenerationIndex;
	TSelectionStrength selectionPressure;
	TKmeansIterations kmeansIterations;
	TNormExponent normExponent;
	TFitnessMode fitnessMode;
	TLogFile logFile;
	Tseed prngSeed;
	TBoolean sortByPosition;
	PyTupleObject*  intergenicCutParamPython;
	TBoolean saveBestIndividual;
	TFitness          nonCodingGenomeFitness;
	TIndexGene genomeSizeLimit;
    TBoolean          shuffleGenome;
	if (PyArg_ParseTuple(args, "OiiiiifffffffiiififisfiOiOOOOfOi",&capsule,&geneSize, &numberOfInitialGenes, &sizeDataBuffer,&sizeDataPoint, &genomeSizeLimit,&largeDuplicationRate,&largeDeletionRate,&largeTranslocationRate,&largeInvertionRate,&pointSubstitutionRate,&pointInsertionRate,&pointDeletionRate,&moreThanOnePointSubstitution,&populationSize,&currentGenerationIndex,&selectionPressure,&kmeansIterations,&normExponent,&fitnessMode,&logFile,&prngSeed,&sortByPosition,&intergenicCutParamPython,&saveBestIndividual ,&mutLawsPython,&initLawsPython,&boundaryConditionsPython,&mutLawsGeneElementPython, &nonCodingGenomeFitness ,&unknownDimentionRandomValuesGenerator, &shuffleGenome)){
		simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
		ParseAllGeneStatisticalLaws(mutLawsPython, simulation->a_parameters->mutationLaws);
		ParseAllGeneStatisticalLaws(initLawsPython,  simulation->a_parameters->initPositionsPDF);
		ParseStatisticalLaw(unknownDimentionRandomValuesGenerator,  simulation->a_parameters->unknownDimentionRandomValuesGenerator);
		ParseAllGeneBoundaryConditions(boundaryConditionsPython,  simulation->a_parameters->boundaryConditions);
		ParseMutationLawGeneElement(mutLawsGeneElementPython,  simulation->a_parameters->geneEleMutationProbsLaw);
		ParseIntergenicCutParam(intergenicCutParamPython, simulation->a_parameters->intergenicCutParameter);
		simulation->a_parameters->geneSize = geneSize;
		simulation->a_parameters->numberOfInitialGenes = numberOfInitialGenes;
		simulation->a_parameters->sizeDataBuffer = sizeDataBuffer;
		simulation->a_parameters->sizeDataPoint = sizeDataPoint;
		simulation->a_parameters->largeDuplicationRate = largeDuplicationRate;
		simulation->a_parameters->largeDeletionRate = largeDeletionRate;
		simulation->a_parameters->largeTranslocationRate = largeTranslocationRate;
		simulation->a_parameters->largeInvertionRate = largeInvertionRate;
		simulation->a_parameters->pointSubstitutionRate = pointSubstitutionRate;
		simulation->a_parameters->pointInsertionRate = pointInsertionRate;
		simulation->a_parameters->pointDeletionRate = pointDeletionRate;
		simulation->a_parameters->moreThanOnePointSubstitution = moreThanOnePointSubstitution;
		simulation->a_parameters->populationSize = populationSize;
		simulation->a_parameters->currentGenerationIndex = currentGenerationIndex;
		simulation->a_parameters->selectionPressure = selectionPressure;
		simulation->a_parameters->kmeansIterations = kmeansIterations;
		simulation->a_parameters->normExponent = normExponent;
		simulation->a_parameters->fitnessMode = fitnessMode;
		simulation->a_parameters->logFile = logFile;
		simulation->a_parameters->prngSeed = prngSeed;
		simulation->a_parameters->sortByPosition = sortByPosition;
		simulation->a_parameters->saveBestIndividual = saveBestIndividual;
		simulation->a_parameters->nonCodingGenomeFitness = nonCodingGenomeFitness;
		simulation->a_parameters->genomeSizeLimit = genomeSizeLimit;
		simulation->a_parameters->shuffleGenome = shuffleGenome;
		return Py_BuildValue("i", 1);
		}
	else{
		return NULL;
		}
	}


static PyObject* SetData(PyObject* self, PyObject* args){
	PyListObject** input;
	CSimulation* simulation;
	PyObject* capsule;
	TBoolean  replace;
	if (!PyArg_ParseTuple(args, "OOi",  &capsule, &input, &replace)) {
		return NULL;
	}
	simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	if (replace == TRUE){
    	simulation->a_startObservationPosInData = 0;
		simulation->a_stopObservationPosInData = 0;
  	}
	TIndexPoint size = (TIndexPoint) PyList_Size((PyObject*)input);
	TIndexCoordinate size_loc = 0;
	TIndexPoint pointToReplace;
	for (TIndexPoint i = 0; i < size; i++) {
		size_loc = (TIndexCoordinate) PyList_Size((PyObject*)PyList_GetItem((PyObject*)input, (Py_ssize_t)i));
		pointToReplace = simulation->SetPoint(size_loc-2);//-1 for label and -1 for cluster
		simulation->a_dataStream->p_arrayPoints->p_elt[pointToReplace].label   = PyInt_AsLong((PyObject*)PyList_GetItem((PyObject*)PyList_GetItem((PyObject*)input, (Py_ssize_t)i), (Py_ssize_t)0));
		simulation->a_dataStream->p_arrayPoints->p_elt[pointToReplace].cluster =  PyInt_AsLong((PyObject*)PyList_GetItem((PyObject*)PyList_GetItem((PyObject*)input, (Py_ssize_t)i), (Py_ssize_t)1));
        for (TIndexCoordinate j = 2; j < size_loc; j++){
			for (TIndexCoordinateElement k = 0;k<COORDINATE_LENGTH;k++){
				simulation->a_dataStream->p_arrayPoints->p_elt[pointToReplace].p_coordinates->p_elt[j-2].coordinate[k] = (TCoordinateElement) PyInt_AsLong((PyObject*)PyList_GetItem((PyObject*)PyList_GetItem((PyObject*)PyList_GetItem((PyObject*)input, (Py_ssize_t)i), (Py_ssize_t)j),(Py_ssize_t)k));//+1for label and +1 for cluster
			}
		}
	}
	return Py_BuildValue("i", 1);
}


static void PieceOfDataToArray(TData* data,PyListObject*** output,TIndexPoint start, TIndexPoint stop, TIndexPoint counter_state){
	TIndexPoint size_loc = 0;
	TIndexPoint counter = counter_state;
	for (TIndexPoint i = start; i < stop; i++) {
		size_loc = data->p_arrayPoints->p_elt[i].p_coordinates->nbElts;
		PyList_SetItem((PyObject*)PyList_GetItem((PyObject*)output, (Py_ssize_t)counter), (Py_ssize_t)0,  PyInt_FromLong(data->p_arrayPoints->p_elt[i].label));
		PyList_SetItem((PyObject*)PyList_GetItem((PyObject*)output, (Py_ssize_t)counter), (Py_ssize_t)1,  PyInt_FromLong(data->p_arrayPoints->p_elt[i].cluster));
		for (TIndexCoordinate j = 2; j < size_loc+2; j++){
			for (TIndexCoordinateElement k = 0;k<COORDINATE_LENGTH;k++){
				PyList_SetItem((PyObject*)PyList_GetItem((PyObject*)PyList_GetItem((PyObject*)output, (Py_ssize_t)counter), (Py_ssize_t)j), (Py_ssize_t)k,   PyFloat_FromDouble(data->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j-2].coordinate[k]));
			}
		}
		counter++;
	}
}

static PyObject* DataToArray(TData* data,TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData, PyListObject*** output_i){
	if (data == NULL){
		return (PyObject*) PyList_New(0);
	}
	PyListObject*** output = output_i;
	TIndexPoint size = data->p_arrayPoints->nbElts;
	if  (startObservationPosInData < stopObservationPosInData) PieceOfDataToArray(data, output, startObservationPosInData, stopObservationPosInData,0);
	else {
		PieceOfDataToArray(data, output,startObservationPosInData, size, 0);
		PieceOfDataToArray(data, output, 0, stopObservationPosInData, size - startObservationPosInData);
	}
	return Py_BuildValue("i", 1);
}

static PyObject* GetData(PyObject* self, PyObject* args){
	PyObject* capsule;
	CSimulation* simulation;
	PyListObject*** output;
	if (!PyArg_ParseTuple(args, "OO", &capsule,&output)){
		return NULL;
		}

	simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	DataToArray(simulation->a_dataStream,simulation->a_startObservationPosInData, simulation->a_stopObservationPosInData, output);
	return Py_BuildValue("i", 1);
	}






static PyObject* LoadData(PyObject* self, PyObject* args){

	PyObject* capsule;
	TBackupFileName backup_file_name;
	if (!PyArg_ParseTuple(args, "sO",&backup_file_name,&capsule)){
		return NULL;
		}
	gzFile * backup_file = (gzFile *) gzopen (backup_file_name, "r");
	if (! backup_file) {
		fprintf (stderr, "gzopen of '%s' failed: %s.\n", backup_file_name,
        strerror (errno));
        exit (EXIT_FAILURE);
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	for (TIndexPoint i = 0;i<simulation->a_dataStream->p_arrayPoints->nbElts;i++){
		gzread((gzFile) backup_file, &simulation->a_dataStream->p_arrayPoints->p_elt[i].label, sizeof(TGeneElement));
		gzread((gzFile) backup_file, &simulation->a_dataStream->p_arrayPoints->p_elt[i].cluster, sizeof(TGeneElement));
		for (TIndexCoordinate j = 0; j<simulation->a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++){
			gzread((gzFile) backup_file, &simulation->a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], sizeof(TCoordinateElement));
			gzread((gzFile) backup_file, &simulation->a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE], sizeof(TCoordinateElement));
			}
		}
	gzclose ((gzFile)backup_file);
	return Py_BuildValue("i", 1);
	}

static PyObject* SaveData(PyObject* self, PyObject* args){

	PyObject* capsule;
	TBackupFileName backup_file_name;
	if (!PyArg_ParseTuple(args, "sO",&backup_file_name,&capsule)){
		return NULL;
		}
	gzFile * backup_file = (gzFile *) gzopen (backup_file_name, "w");
	if (! backup_file) {
		fprintf (stderr, "gzopen of '%s' failed: %s.\n", backup_file_name,
        strerror (errno));
        exit (EXIT_FAILURE);
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	for (TIndexPoint i = 0;i<simulation->a_dataStream->p_arrayPoints->nbElts;i++){
		gzwrite((gzFile) backup_file, &simulation->a_dataStream->p_arrayPoints->p_elt[i].label, sizeof(TGeneElement));
		gzwrite((gzFile) backup_file, &simulation->a_dataStream->p_arrayPoints->p_elt[i].cluster, sizeof(TGeneElement));
		for (TIndexCoordinate j = 0; j<simulation->a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->nbElts; j++){
			gzwrite((gzFile) backup_file, &simulation->a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[DIMENSION_IN_COORDINATE], sizeof(TCoordinateElement));
			gzwrite( (gzFile)backup_file, &simulation->a_dataStream->p_arrayPoints->p_elt[i].p_coordinates->p_elt[j].coordinate[POSITION_IN_COORDINATE], sizeof(TCoordinateElement));
			}
		}
	gzclose ((gzFile)backup_file);
	return Py_BuildValue("i", 1);
	}

static PyObject* IterateASimulation_traductor(PyObject* self, PyObject* args){
	CSimulation* simulation;
	PyObject* capsule;
	TGenerationsIndex nbIterations;
	if (!PyArg_ParseTuple(args, "Oi", &capsule,& nbIterations)){
		return NULL;
		}
	simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	//simulation->PrintGenomes();
	simulation->EvolveDuringTsteps(nbIterations);
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* ResetASimulationPopulation_traductor(PyObject* self, PyObject* args){
	CSimulation* simulation = SimulationPythonToC(args);
	simulation->InitRandomPopulation();
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* DeleteSimulation_traductor(PyObject* self, PyObject* args){
	CSimulation*  simulation = SimulationPythonToC(args);
	delete simulation;
	Py_INCREF(Py_None);
	return Py_None;
	}
/*
static PyListObject* IndividualStatistics2PyList(TIndividualFeatures individualFeatures){
		PyListObject* list =  (PyObject*) PyList_New(0);
		PyList_Append((PyObject*)list, PyFloat_FromDouble(individualFeatures.fitness));
		PyList_Append((PyObject*)list, PyFloat_FromDouble(individualFeatures.codingRatio));
		PyList_Append((PyObject*)list, PyFloat_FromDouble(1.*individualFeatures.genomeSize));
		return list;
	}
*/

static PyObject* Printgenomes(PyObject* self, PyObject* args){
	CSimulation*  simulation = SimulationPythonToC(args);
	simulation->PrintGenomes();
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* ClassifyData(PyObject*, PyObject* args){
	CSimulation* simulation;
	PyObject* capsule;
	TIndexIndividual individualIndex;
	PyListObject*** output;
	if (!PyArg_ParseTuple(args, "OiO", &capsule,&individualIndex,&output)){
		return NULL;
		}
	simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	CIndividual* individual = simulation->a_simulationState->r_population->a_populationPresent.p_elt[individualIndex];
	switch(simulation->a_parameters->fitnessMode){
		case FITNESS_MODE_CLUSTERING:{
			CFitnessCluster* fitnessClassifier = (CFitnessCluster*) simulation->a_fitnessEvaluator;
			DataToArray(((CFitnessCluster*)fitnessClassifier)->ClassifyData( individual->a_genome,individual->a_geneSize, simulation->a_dataStream,simulation->a_startObservationPosInData, simulation->a_stopObservationPosInData),simulation->a_startObservationPosInData, simulation->a_stopObservationPosInData,output);
			return Py_BuildValue("i", 1);
		}
		case FITNESS_MODE_CLUSTERING_KMEANS:{
			CFitnessClusterKmeans* fitnessClassifier = (CFitnessClusterKmeans*) simulation->a_fitnessEvaluator;
			DataToArray(((CFitnessClusterKmeans*)fitnessClassifier)->ClassifyData( individual->a_genome,individual->a_geneSize, simulation->a_dataStream,simulation->a_startObservationPosInData, simulation->a_stopObservationPosInData),simulation->a_startObservationPosInData, simulation->a_stopObservationPosInData,output);
			return Py_BuildValue("i", 1);
		}
	}
	Py_INCREF(Py_None);
	return Py_None;
}


void IndividualCapsuleDestructor(PyObject* capsule){
	CIndividual* individual = (CIndividual*) PyCapsule_GetPointer(capsule,NAME_CAPSULE_INDIVIDUAL);
	delete individual;
	}

static PyObject* GetIndividual_traductor(PyObject* self, PyObject* args){
	TIndexIndividual individualIndex;
	CSimulation* simulation;
	PyObject* capsule;
	if (!PyArg_ParseTuple(args, "Oi", &capsule,&individualIndex)){
		return NULL;
		}
	simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	CIndividual* individual_to_be_copied   = simulation->a_simulationState->r_population->a_populationPresent.p_elt[individualIndex];
	CIndividual* individual_to_be_returned = new CIndividual(individual_to_be_copied->a_genome->p_arrayGenes->nbElts, individual_to_be_copied->a_genome->p_arrayGenes->maxNbElts, individual_to_be_copied->a_geneSize, individual_to_be_copied->a_position, individual_to_be_copied->a_p_prng);
	individual_to_be_returned->CopyOtherIndividual(individual_to_be_copied);
	PyObject* capsule_individual = PyCapsule_New(individual_to_be_returned,NAME_CAPSULE_INDIVIDUAL, IndividualCapsuleDestructor);
	return capsule_individual;
	}

static PyObject* SetIndividual_traductor(PyObject* self, PyObject* args){
	TIndexIndividual individualIndex;
	CSimulation* simulation;
	PyObject* capsule;
	PyObject* capsule_individual;
	if (!PyArg_ParseTuple(args, "OOi", &capsule,&capsule_individual,&individualIndex)){
		return NULL;
	}
	simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	CIndividual* individual_to_be_copied   = (CIndividual*) PyCapsule_GetPointer(capsule_individual,NAME_CAPSULE_INDIVIDUAL);
	CIndividual* individual_to_be_replaced = simulation->a_simulationState->r_population->a_populationPresent.p_elt[individualIndex];
	individual_to_be_replaced->CopyOtherIndividual(individual_to_be_copied);
	Py_INCREF(Py_None);
	return Py_None;
}

static void SetOrganismGenotype(PyListObject** genomeInput, CIndividual* individual ){
	TIndexGene size = (TIndexGene) PyList_Size((PyObject*)genomeInput);
	individual->a_genome->p_arrayGenes->nbElts = size;
	for (TIndexGene i = 0; i < size; i++) {
		TIndexGeneElement size_loc = (TIndexGeneElement) PyList_Size((PyObject*)PyList_GetItem((PyObject*)genomeInput, (Py_ssize_t)i));
		if (size_loc < GENE_SIZE){
			for (TIndexGeneElement j = 0; j < size_loc; j++){
				TGeneElement value = (TGeneElement) PyInt_AsLong((PyObject*)PyList_GetItem((PyObject*)PyList_GetItem((PyObject*)genomeInput, (Py_ssize_t)i),(Py_ssize_t)j));
				individual->SetGeneElement(i, j, value);
			}
		}
	}
}

static PyObject* SetIndividualGenotype(PyObject*, PyObject* args){
	PyListObject** genomeInput = (PyListObject**)PyList_New(0);
	PyObject* capsule;
	TIndexIndividual individualIndex;
	if (!PyArg_ParseTuple(args, "OiO", &capsule,&individualIndex,&genomeInput)){
		return NULL;
	}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	CIndividual* individual = simulation->a_simulationState->r_population->a_populationPresent.p_elt[individualIndex];
	SetOrganismGenotype(genomeInput, individual );
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* ChangeSeed(PyObject*, PyObject* args){
	PyObject* capsule;
	Tseed  newSeed;
	if (!PyArg_ParseTuple(args, "Of", &capsule,&newSeed)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	//printf("newsedd %lg\n",newSeed);
	simulation->a_parameters->prngSeed = newSeed;
	simulation->ChangeSeed(newSeed);
	Py_INCREF(Py_None);
	return Py_None;
	}
static PyObject*  Resetdata (PyObject*, PyObject* args){
	PyObject* capsule;
	TIndexPoint  memoryLength;
	if (!PyArg_ParseTuple(args, "Oi", &capsule,&memoryLength)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	//printf("newsedd %lg\n",newSeed);
	simulation->a_startObservationPosInData = 0;
	simulation->a_stopObservationPosInData  = 0;
	if (simulation->a_dataStream->p_arrayPoints->nbElts<memoryLength){
		printf("error new data memory length bigger than older one!!!!!!!!!!!!!\n");
		}
	simulation->a_dataStream->p_arrayPoints->nbElts = memoryLength;
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* SaveClassifyDataSetIndiv(PyObject*, PyObject* args){
	CSimulation* simulation;
	PyObject* capsule;
	TIndexIndividual individualIndex;
	TFileName  nameFile;
	if (!PyArg_ParseTuple(args, "Ois", &capsule,&individualIndex,&nameFile)){
		return NULL;
		}
	simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	simulation->SaveIndividualClassifiedData(simulation->a_dataStream, simulation->a_startObservationPosInData, simulation->a_stopObservationPosInData, nameFile, individualIndex);
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* SaveStats(PyObject*, PyObject* args){
	PyObject* capsule;
	TFileName  nameFile;
	if (!PyArg_ParseTuple(args, "Os", &capsule,&nameFile)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	simulation->SavePopulationFeatures(nameFile);
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* SaveGenotype(PyObject*, PyObject* args){
	PyObject* capsule;
	TIndexIndividual  individual;
	TFileName  nameFile;
	if (!PyArg_ParseTuple(args, "Ois", &capsule,&individual,&nameFile)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	simulation->SaveIndividualGenome(nameFile,individual);
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* SavePhenotype(PyObject*, PyObject* args){
	PyObject* capsule;
	TIndexIndividual  individual;
	TFileName  nameFile;
	if (!PyArg_ParseTuple(args, "Ois", &capsule,&individual,&nameFile)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	simulation->SaveIndividualPhenotype(nameFile,individual);
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* SaveIndividualStats(PyObject*, PyObject* args){
	PyObject* capsule;
	TIndexIndividual  individual;
	TFileName  nameFile;
	if (!PyArg_ParseTuple(args, "Ois", &capsule,&individual,&nameFile)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	simulation->SaveIndividualFeatures(nameFile,  individual);
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* ModifyPseudogenesAllPopulation(PyObject*, PyObject* args){
	PyObject* capsule;
	PyTupleObject** modificationLawsPython = (PyTupleObject**)PyTuple_New(0);

	if (!PyArg_ParseTuple(args, "OO", &capsule,&modificationLawsPython)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	TMutationLaw ** modificationDistributionLawPoints = (TMutationLaw **)  malloc( simulation->a_parameters->geneSize *sizeof(TMutationLaw*));
	for (TIndexGeneElement i = 0; i<simulation->a_parameters->geneSize;i++){
		modificationDistributionLawPoints[i] = (TMutationLaw *)  malloc(sizeof(TMutationLaw));
		modificationDistributionLawPoints[i]->mutationLaw = (TMutationLawParameters*)  malloc(sizeof(TMutationLawParameters));
		modificationDistributionLawPoints[i]->mutationLaw->p_transitionMatrix = (TTransitionMatrix*) malloc(sizeof(TTransitionMatrix));
		}//the free is done in the simulation
	ParseAllGeneStatisticalLaws(modificationLawsPython, modificationDistributionLawPoints);
	simulation->ModifyPseudogenesAllPopulation(modificationDistributionLawPoints);
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* ModifyPseudogenesOneIndividual(PyObject*, PyObject* args){
	PyObject* capsule;
	PyTupleObject** modificationLawsPython = (PyTupleObject**)PyTuple_New(0);
	TIndexIndividual  individual;
	if (!PyArg_ParseTuple(args, "OOi", &capsule,&modificationLawsPython,&individual)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	TMutationLaw ** modificationDistributionLawPoints = (TMutationLaw **)  malloc( simulation->a_parameters->geneSize *sizeof(TMutationLaw*));
	for (TIndexGeneElement i = 0; i<simulation->a_parameters->geneSize;i++){
		modificationDistributionLawPoints[i] = (TMutationLaw *)  malloc(sizeof(TMutationLaw));
		modificationDistributionLawPoints[i]->mutationLaw = (TMutationLawParameters*)  malloc(sizeof(TMutationLawParameters));
		modificationDistributionLawPoints[i]->mutationLaw->p_transitionMatrix = (TTransitionMatrix*) malloc(sizeof(TTransitionMatrix));
		}//the free is done in the simulation
	ParseAllGeneStatisticalLaws(modificationLawsPython, modificationDistributionLawPoints);
	simulation->ModifyPseudogenesOneIndividual(modificationDistributionLawPoints,individual);
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* OrderGenomeOneIndividual(PyObject*, PyObject* args){
	PyObject* capsule;
	TIndexIndividual  individual;
	if (!PyArg_ParseTuple(args, "Oi", &capsule,&individual)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	simulation->OrderGenomeOneIndividual(individual);
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* OrderGenomeAllPopulation(PyObject*, PyObject* args){
	CSimulation*  simulation = SimulationPythonToC(args);
	simulation->OrderGenomesAllPopulation();
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* ShuffleGenomeAllPopulation(PyObject*, PyObject* args){
	CSimulation*  simulation = SimulationPythonToC(args);
	simulation->ShuffleGenomeAllPopulation();
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* ShuffleGenomeOneIndividual(PyObject*, PyObject* args){
	PyObject* capsule;
	TIndexIndividual  individual;
	if (!PyArg_ParseTuple(args, "Oi", &capsule,&individual)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	simulation->ShuffleGenomeOneIndividual( individual);
	Py_INCREF(Py_None);
	return Py_None;
	}
	
static 	PyObject* getGenomeSize(PyObject*, PyObject* args){
	PyObject* capsule;
	TIndexIndividual  individualIndex;
	if (!PyArg_ParseTuple(args, "Oi", &capsule,&individualIndex)){
		return NULL;
		}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	CIndividual* individual = simulation->a_simulationState->r_population->a_populationPresent.p_elt[individualIndex];
	return Py_BuildValue("i", individual->a_genome->p_arrayGenes->nbElts);
	}

static PyObject* getIndividualFeaturesSize(PyObject*, PyObject* args){
	PyObject* capsule;
	if (!PyArg_ParseTuple(args, "O", &capsule)){
		return NULL;
	}
	return Py_BuildValue("i", NB_INDIVIDUAL_FEATURES);
}

static PyObject* ReturnIndividualGenome(PyObject* self, PyObject* args){
	PyListObject** output;
	CSimulation* simulation;
	PyObject* capsule;
	TIndexIndividual individualIndex;
	if (!PyArg_ParseTuple(args, "OiO", &capsule, &individualIndex,&output)){
		return NULL;
		}
	simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	CIndividual* individual = simulation->a_simulationState->r_population->a_populationPresent.p_elt[individualIndex];
	for (TIndexGene i = 0; i < individual->a_genome->p_arrayGenes->nbElts; i++) {
		for (TIndexGeneElement j = 0; j < individual->a_geneSize; j++){
			PyList_SetItem(PyList_GetItem((PyObject*)output, (Py_ssize_t)i),(Py_ssize_t)j,  PyFloat_FromDouble(1.*individual->a_genome->p_arrayGenes->p_elt[i].gene[j]));
		}
	}
	Py_INCREF(Py_None);
	return Py_None;
	}

static PyObject* ReturnIndividualStatistics(PyObject* self, PyObject* args){
	PyObject* capsule;
	PyListObject** output;
	if (!PyArg_ParseTuple(args, "OO", &capsule,&output)){
		return NULL;
	}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	TIndexIndividual size_population  = simulation->a_simulationState->r_population->a_populationSize;
	TIndexIndividual size_python_list = (TIndexIndividual) PyObject_Length( (PyObject*) output);
	if (size_population != size_python_list){
		return Py_BuildValue("i", 0);
    }
  	if (PyObject_Length( (PyObject*) PyList_GetItem((PyObject*)output, (Py_ssize_t)0)) < NB_INDIVIDUAL_FEATURES){
    	return Py_BuildValue("i", 0);
    }
	for (TIndexIndividual i = 0; i < size_python_list; i++) {
		TIndividualFeatures individualFeatures = simulation->a_simulationState->r_population->a_populationPresent.p_elt[i]->getStatistics();
		PyList_SetItem(PyList_GetItem((PyObject*)output, (Py_ssize_t)i),(Py_ssize_t)FEATURE_FITNESS_POS,  PyFloat_FromDouble(individualFeatures.fitness*1.0));
		PyList_SetItem(PyList_GetItem((PyObject*)output, (Py_ssize_t)i),(Py_ssize_t)FEATURE_CODING_RATIO_POS,  PyFloat_FromDouble(individualFeatures.codingRatio*1.0));
		PyList_SetItem(PyList_GetItem((PyObject*)output, (Py_ssize_t)i),(Py_ssize_t)FEATURE_GENOME_LENGTH_POS,  PyFloat_FromDouble(individualFeatures.genomeSize*1.0));
		PyList_SetItem(PyList_GetItem((PyObject*)output, (Py_ssize_t)i),(Py_ssize_t)FEATURE_FITNESS_INDEX_POS,  PyFloat_FromDouble(individualFeatures.fitnessIndex*1.0));
		PyList_SetItem(PyList_GetItem((PyObject*)output, (Py_ssize_t)i),(Py_ssize_t)FEATURE_PARENT_INDEX_POS,  PyFloat_FromDouble(individualFeatures.parentIndex*1.0));
		PyList_SetItem(PyList_GetItem((PyObject*)output, (Py_ssize_t)i),(Py_ssize_t)FEATURE_INDIV_INDEX_POS,  PyFloat_FromDouble(individualFeatures.individualIndex*1.0));
   		//printf("%lg %lg %i %i %i %i\n",individualFeatures.fitness,individualFeatures.codingRatio,individualFeatures.genomeSize,individualFeatures.parentIndex,individualFeatures.individualIndex);
		}
	Py_INCREF(Py_None);
	return Py_None;
	}

static void ComputeLenLinesPieceOfDataToArray(TData* data,PyListObject* output,TIndexPoint start, TIndexPoint stop, TIndexPoint counter_state){
	TIndexPoint size_loc = 0;
	TIndexPoint counter = counter_state;
	for (TIndexPoint i = start; i < stop; i++) {
		size_loc = data->p_arrayPoints->p_elt[i].p_coordinates->nbElts;
		PyList_SetItem((PyObject*)output, (Py_ssize_t)counter,  PyInt_FromLong(size_loc + 2));
		counter++;
	}
}

static PyObject* ComputeLenLinesData(TData* data,TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData, PyListObject* lenLines_i){
	if (data == NULL){
		return Py_BuildValue("i", 0);
	}
	PyListObject* output = lenLines_i;
	TIndexPoint size = data->p_arrayPoints->nbElts;
	if  (startObservationPosInData < stopObservationPosInData){
		ComputeLenLinesPieceOfDataToArray(data, output, startObservationPosInData, stopObservationPosInData,0);
	}
	else {
		ComputeLenLinesPieceOfDataToArray(data, output,startObservationPosInData, size, 0);
		ComputeLenLinesPieceOfDataToArray(data, output, 0, stopObservationPosInData, size - startObservationPosInData);
	}
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* GetIndividualPhenotype(PyObject*, PyObject* args){
	CSimulation* simulation;
	PyObject* capsule;
	PyListObject*** output;
	PyListObject* lenLines;
	TIndexIndividual individualIndex;

	if (!PyArg_ParseTuple(args, "OiOO", &capsule,&individualIndex,&output,&lenLines)){
		return NULL;
		}
	simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	CIndividual* individual = simulation->a_simulationState->r_population->a_populationPresent.p_elt[individualIndex];
	TData* phenotype = simulation->a_fitnessEvaluator->ComputePhenotype( individual->a_genome, individual->a_geneSize);
	DataToArray(phenotype,0, phenotype->p_arrayPoints->nbElts,output);
	ComputeLenLinesData(phenotype,0, phenotype->p_arrayPoints->nbElts,lenLines);
	return Py_BuildValue("i", phenotype->p_arrayPoints->nbElts);
	}

static PyObject* ComputeFitnesses(PyObject*, PyObject* args){
	PyObject* capsule;
	if (!PyArg_ParseTuple(args, "O", &capsule)){
		return NULL;
	}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	simulation->a_simulationState->r_population->ComputeIndividualFitness(simulation->a_dataStream, simulation->a_startObservationPosInData, simulation->a_stopObservationPosInData, (CFitnessEvaluator*) simulation->a_fitnessEvaluator,simulation->a_parameters->saveBestIndividual);
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* ModifyOneGene(PyObject*,PyObject* args){
    PyObject* capsule;
    TIndexIndividual individualIndex;
    TIndexGene geneIndex;
    PyListObject* newGene;// = PyTuple_GetItem((PyObject*)geneEleMutationProbsLawPython, (Py_ssize_t)0);
    if (!PyArg_ParseTuple(args, "OiiO",&capsule,&individualIndex,&geneIndex,&newGene)){
        return NULL;
    }
    CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
    CIndividual* indiv = simulation->a_simulationState->r_population->a_populationPresent.p_elt[individualIndex];
    for(TGeneElement j = 0; j < indiv->a_geneSize; j++){
        indiv->a_genome->p_arrayGenes->p_elt[geneIndex].gene[j] = PyInt_AsLong((PyObject*)PyList_GetItem((PyObject*)newGene,(Py_ssize_t)j));
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* ComputeIndividualFitness(PyObject*, PyObject* args){
	PyObject* capsule;
	TIndexIndividual  individualIndex;
	if (!PyArg_ParseTuple(args, "Oi", &capsule,&individualIndex)){
		return NULL;
	}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	simulation->a_simulationState->r_population->a_populationPresent.p_elt[individualIndex]->ComputeFitness(simulation->a_dataStream, simulation->a_startObservationPosInData, simulation->a_stopObservationPosInData, (CFitnessEvaluator*) simulation->a_fitnessEvaluator);
	Py_INCREF(Py_None);
	return Py_BuildValue("f", simulation->a_simulationState->r_population->a_populationPresent.p_elt[individualIndex]->a_fitness);
}

static PyObject* GetMaximalGeneSize(PyObject*, PyObject* args){
	return Py_BuildValue("i", GENE_SIZE);
	}

static PyObject* SetBestIndividualGenotype(PyObject*, PyObject* args){
	PyObject* capsule;
	PyListObject** genomeInput = (PyListObject**)PyList_New(0);
	if (!PyArg_ParseTuple(args, "OO", &capsule,&genomeInput)){
		return NULL;
	}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	CIndividual* individual =  simulation->a_simulationState->r_population->a_bestOldIdividual;
	SetOrganismGenotype(genomeInput, individual );
	Py_INCREF(Py_None);
	return Py_None;
}


static PyObject* SetSizeDataBuffer(PyObject*, PyObject* args){
    PyObject* capsule;
	TIndexPoint  sizeDataBuffer;
	TIndexCoordinate sizeDataPoint;
	if (!PyArg_ParseTuple(args, "iiO", &sizeDataBuffer,&sizeDataPoint,&capsule)){
		return NULL;
	}
	CSimulation* simulation = (CSimulation*) PyCapsule_GetPointer(capsule,NAME_CAPSULE);
	simulation->a_parameters->sizeDataBuffer	= sizeDataBuffer;
	simulation->a_parameters->sizeDataPoint	= sizeDataPoint;
	simulation->ReallocateDataStreamMem();
	Py_INCREF(Py_None);
	return Py_None;
  }
  
static PyMethodDef modevoevo_funcs[] = {
    {"esimulation",(PyCFunction)CreatSimulation_traductor,METH_VARARGS, NULL},
    {"eiterate", IterateASimulation_traductor, METH_VARARGS, NULL},
    {"ereset", ResetASimulationPopulation_traductor, METH_VARARGS, NULL},
    {"edelete",DeleteSimulation_traductor, METH_VARARGS, NULL},
    {"esetdata",SetData,METH_VARARGS, NULL},
    {"egetstats",ReturnIndividualStatistics,METH_VARARGS, NULL},
    {"egetindividualgenome",ReturnIndividualGenome,METH_VARARGS, NULL},
    {"egetindividualphenotype",GetIndividualPhenotype,METH_VARARGS, NULL},
    {"esetindividualgenome",SetIndividualGenotype,METH_VARARGS, NULL},
    {"esetbestindividualgenome",SetBestIndividualGenotype,METH_VARARGS, NULL},
    {"eloaddata",LoadData,METH_VARARGS, NULL},
    {"esavedata",SaveData,METH_VARARGS, NULL},
    {"esavesimulation",SaveSimulation,METH_VARARGS, NULL},
    {"eloadsimulation",InitSimulationFromBackup,METH_VARARGS, NULL},
    {"eclassifydata",ClassifyData,METH_VARARGS, NULL},
    {"eprintgenomes",Printgenomes,METH_VARARGS, NULL},
    {"egetparameters",Getparameters,METH_VARARGS, NULL},
    {"echangeseed",ChangeSeed,METH_VARARGS, NULL},
    {"eresetdata",Resetdata,METH_VARARGS, NULL},
    {"esavestats",SaveStats,METH_VARARGS, NULL},
	{"esavegenotype",SaveGenotype,METH_VARARGS, NULL},
	{"esavephenotype",SavePhenotype,METH_VARARGS, NULL},
	{"esaveindividualstat",SaveIndividualStats,METH_VARARGS, NULL},
	{"eclassifydataset",SaveClassifyDataSetIndiv,METH_VARARGS, NULL},
	{"emodifypseudogenesallpopulation",ModifyPseudogenesAllPopulation,METH_VARARGS,NULL},
	{"emodifypseudogenesoneindividual",ModifyPseudogenesOneIndividual,METH_VARARGS,NULL},
	{"eordergenomeindividual",OrderGenomeOneIndividual,METH_VARARGS,NULL},
	{"eordergenomeallpopulation",OrderGenomeAllPopulation,METH_VARARGS,NULL},
	{"eshuffleallpopulationgenome",ShuffleGenomeAllPopulation,METH_VARARGS,NULL},
	{"eshuffleoneindividualgenome",ShuffleGenomeOneIndividual,METH_VARARGS,NULL},
	{"egetgenomesize",getGenomeSize,METH_VARARGS,NULL},
	{"egetindividualfeaturesnb",getIndividualFeaturesSize,METH_VARARGS,NULL},
	{"egetdata",GetData,METH_VARARGS,NULL},
    {"ecomputefitnesses",ComputeFitnesses,METH_VARARGS,NULL},
    {"ecomputeindividualfitness",ComputeIndividualFitness,METH_VARARGS,NULL},
    {"egetmaximalgenesize",GetMaximalGeneSize,METH_VARARGS,NULL},
    {"egetindividual",GetIndividual_traductor,METH_VARARGS,NULL},
    {"esetindividual",SetIndividual_traductor,METH_VARARGS,NULL},
    {"esetsizedatabuffer",SetSizeDataBuffer,METH_VARARGS,NULL},
    {"esetparameters",SetParameters,METH_VARARGS,NULL},
    {"emodifygene",ModifyOneGene,METH_VARARGS,NULL},
    {NULL, NULL, 0, NULL}
};
void initmodevoevo_c(void){
    Py_InitModule3("modevoevo_c", modevoevo_funcs,
                   "modevoevo module");
}



#ifdef __cplusplus
}  // extern "C"
#endif
