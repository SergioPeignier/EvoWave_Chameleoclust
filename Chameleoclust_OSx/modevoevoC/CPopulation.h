#include <zlib.h>
#include <errno.h>
#include <string.h>
#include "evoevo_def_state_pop.h"
#include "evoevo_def_base.h"
#include "evoevo_def_parameters.h"
#include "CPrng.h"
#include "CFitnessEvaluator.h"
#ifndef POP
#define POP
class CPopulation{
	
	public :
	CIndividual*        a_bestOldIdividual;
	TArrayIndividuals   a_populationPresent;
	TReproductionProb*  a_p_reproductionProbs;
	TPopulationSize     a_populationSize;
	TSelectionStrength  a_selectionPressure;
	TIndexIndividual*   perIndividualOffspring;
	TIndexIndividual*   a_availablePositions;
	TIndexIndividual*   a_occupiedPositions;
	TIndexIndividual    a_nbAvailablePositions;
	TIndexIndividual    a_nbOccupiedPositions;
	TGenerationsIndex   a_paramCurentGenerationindex;
	CPrng*              a_p_prng;
	
	CPopulation();
	CPopulation(TPopulationSize populationSize, TSelectionStrength selectionPressure, CPrng* prng, TGenerationsIndex paramCurentGenerationIndex);
	~CPopulation();
	CPopulation(gzFile* backup_file , CPrng* prng);//
	void loadIndividualsArray(gzFile* backup_file,CPrng* prng);//
	void save( gzFile* backup_file );//
	void saveIndividualsArray(gzFile* backup_file);//
	TIndexIndividual* SelectedFuturProgeny();
	void GenerateReproductionProbas();
	TMutationRate ReproductionProba(TIndexIndividual rank);
	void GiveBirth (TIndexIndividual* fathersIndexes);
	void GenerateRandomPopulation(TIndexGene initNumberGenes,TIndexGene maxNbGenes,TIndexGeneElement geneSize,TMutationLaw** initialPositionsPDF, TBoundaryConditions** geneBoundaryConditions );
	void GenerateRandomPopulation(TIndexGene initNumberGenes,TIndexGene maxNbGenes,TIndexGeneElement geneSize,TMutationLaw** initialPositionsPDF , TBoundaryConditions** geneBoundaryConditions,  TIndexIndividual initGene, TIndexIndividual endGene);
	static int CompareIndividualsByFitness(const void* A_i, const void* B_i);
	static int CompareIndividualsByPosition(const void* A_i, const void* B_i);
	void PrintAllGenomes();
	void checkAvailableOccupiedPositions(TIndexIndividual* perIndividualOffspring);
	void MutateNextGeneration(TSimulationParameters * simulationParameters);
	void SortNewGeneration();
	void ComputeNextGeneration(TSimulationParameters * simulationParameters);
	void ComputeIndividualFitness( TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData,  CFitnessEvaluator * fitnessComputer,TBoolean saveBestIndividual);
	void EvolveDuringTsteps(TSimulationParameters * simulationParameters,TGenerationsIndex T, TData* dataStream,TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData, CFitnessEvaluator * fitnessComputer);
	void SortByFitness();
	void ShuffleGenomes();
	void SetPRNG(CPrng* new_prng);
	
	void SaveIndividualGenome(TFileName nameSaveFile,TIndexIndividual individual);
	void SaveIndividualFeatures(TFileName nameSaveFile, TIndexIndividual individual);
	void SaveIndividualPhenotype( CFitnessEvaluator * fitnessComputer, TFileName nameSaveFile,TIndexIndividual individual);
	void SavePopulationFeatures(TFileName nameSaveFile);
	void SaveIndividualClassifiedData(TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData ,CFitnessEvaluator * fitnessComputer, TFileName nameSaveFile,TIndexIndividual individual);

	void ModifyPseudogenesAllPopulation(TMutationLaw** newPseudogenePDF, TBoundaryConditions** geneBoundaryConditions);
	void ModifyPseudogenesOneIndividual(TMutationLaw** newPseudogenePDF, TBoundaryConditions** geneBoundaryConditions,TIndexIndividual individual);
	
	void OrderGenomesAllPopulation();
	void OrderGenomeOneIndividual(TIndexIndividual individual);
	
	void ShuffleGenomeAllPopulation();
	void ShuffleGenomeOneIndividual(TIndexIndividual individual);
	
	TIndexIndividual GetBestIndividualIndex();
	void SaveOrUpdateBestIndividual(TBoolean saveBestIndividual);

};
#endif
