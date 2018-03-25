#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <zlib.h>
#include <errno.h>
#include <string.h>
#include "evoevo_def_state_pop.h"
#include "evoevo_def_base.h"
#include "evoevo_def_parameters.h"
#include "CPrng.h"
#include "CFitnessModeTest1.h"
#include "CFitnessRandom.h"
#include "CFitnessCluster.h"
#include "CFitnessClusterKmeans.h"
#include "CFitnessEvaluator.h"

#ifndef SIM
#define SIM
class CSimulation{
	public :
	TSimulationState* a_simulationState;
	CPrng* a_p_prng;
	TSimulationParameters* a_parameters;
	TData* a_dataStream;
	CFitnessEvaluator* a_fitnessEvaluator; //changer avec un heritage
	TIndexPoint a_startObservationPosInData;
	TIndexPoint a_stopObservationPosInData;
	CSimulation(TSimulationParameters*);
	CSimulation(gzFile* backup_file,FILE* backup_state_rng, TSimulationParameters* parameters);//
	void save( gzFile* backup_file , FILE* backup_state_rng);//
	~CSimulation();
	void InitRandomPopulation();
	void ChangePopulationSize(TPopulationSize newPopulationSize);
	void ChangeSelectionPressure( TSelectionStrength newSelectionPressure);
	TGenome* GetGenomeBestIndivLastGeneration();
	void LoadParameters(TSimulationParameters* parameters);
	void PrintGenomes();
	void EvolveDuringTsteps(TGenerationsIndex T);
	void AllocateDataStreamMem();
	TIndexPoint SetPoint(TIndexCoordinate pointDimentionality);
	void SetPointCoordinate(TCoordinateElement* coordinate, TIndexCoordinate coordinateIndexInPoint, TIndexPoint pointIndex );
	void  PrintDataSet();
	void InitializeFitnessFunction();
	void freeMutationLawsArray(TMutationLaw** mutationLaws);
	void freeBoundaryConditions(TBoundaryConditions** boundaryConditions);
	void freeGeneEleMutationProbLaw(TGeneEleMutationProbsLaw* geneEleMutationProbsLaw);
	void freeDataStreamMem(TData* dataStream);
	void freeMutationLaw(TMutationLaw* mutationLaw);
	void ChangeSeed(Tseed  newSeed);

	void SaveIndividualGenome(TFileName nameSaveFile,TIndexIndividual individual);
	void SaveIndividualFeatures(TFileName nameSaveFile, TIndexIndividual individual);
	void SaveIndividualPhenotype( TFileName nameSaveFile,TIndexIndividual individual);
	void SavePopulationFeatures(TFileName nameSaveFile);
	void SaveIndividualClassifiedData(TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData , TFileName nameSaveFile,TIndexIndividual individual);


	void ModifyPseudogenesAllPopulation(TMutationLaw** newPseudogenePDF);
	void ModifyPseudogenesOneIndividual(TMutationLaw** newPseudogenePDF,TIndexIndividual individual);

	void OrderGenomesAllPopulation();
	void OrderGenomeOneIndividual(TIndexIndividual individual);

	void ShuffleGenomeAllPopulation();
	void ShuffleGenomeOneIndividual(TIndexIndividual individual);
	void ReallocateDataStreamMem();
	void freeIntergenicCutParameter(TIntergenicCut* intergenicCutParameter);

};
#endif
