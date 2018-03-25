#include <zlib.h>
#include <errno.h>
#include <string.h>
#include "evoevo_def_state_pop.h"
#include "evoevo_def_base.h"
#include "evoevo_def_parameters.h"
#include "CPrng.h"
#include "CFitnessEvaluator.h"

#ifndef INDIV
#define INDIV
class CIndividual{
	public :
		TIndexIndividual   a_position;
		TIndexIndividual   a_parentPosition = 0;
		TIndexIndividual   a_fitnessIndex   = 0; //may be erased latter
		TCodingRatio       a_codingRatio;
		TFitness           a_fitness;
		TGenome*           a_genome;
		CPrng*             a_p_prng;
		TIndexGeneElement  a_geneSize;
		TGene * a_p_tempGen;
		TIndexGene nbLargeMutations [NB_LARGE_MUTATIONS];
		
		CIndividual();
		CIndividual(TIndexGene initNumGenes,TIndexGene maxNbGenes, TIndexGeneElement geneSize, TIndexIndividual position, CPrng* prng);
		~CIndividual();
		
		TIndividualFeatures getStatistics();
		TIndexGene LargeCut(TIndexGene cutGeneStart, TIndexGene cutGeneEnd);
		TGeneElement GeneElementRandomMutationWalk(TMutationLaw* mutationPDF, TIndexGeneElement index);
		TGeneElement ComputeGeneElementAfterSubstitutionWithBoundary(TGeneElement oldGeneElement, TGeneElement randomWalk, TBoundaryConditions* geneBoundaryConditions);
		
		void CSetIndividualGenome(TGenome* other);

		CIndividual(gzFile* backup_file , CPrng* prng);//
		void loadGenome(gzFile* backup_file);//
		void save( gzFile* backup_file );//
		void saveGenome(gzFile* backup_file);//
		
		void InitRandom(TMutationLaw** initialPositionsPDF ,TBoundaryConditions** geneBoundaryConditions );
		void LargeDuplication(TIndexGene cutGeneStart, TIndexGene cutGeneEnd, TIndexGene insertionGene, TIntergenicCut* intergenicCutParameter);
		void LargeDuplication(TIntergenicCut* intergenicCutParameter);
		void LargeDeletion(TIndexGene cutGeneStart, TIndexGene cutGeneEnd,TIntergenicCut* intergenicCutParameter);
		void LargeDeletion(TIntergenicCut* intergenicCutParameter);
		void LargeTranslocation(TIndexGene cutGeneStart, TIndexGene cutGeneEnd, TIndexGene insertionGene, TIntergenicCut* intergenicCutParameter);
		void LargeTranslocation(TIntergenicCut* intergenicCutParameter);
		void LargeInvertion(TIndexGene cutGeneStart, TIndexGene cutGeneEnd,TIntergenicCut* intergenicCutParameter);
		void LargeInvertion(TIntergenicCut* intergenicCutParameter);
		void PrintGenome();
		void DoAllLargeMutations(TSimulationParameters * simulationParameters);
		void DoGenePointSubstitution(TMutationLaw** mutationPDF,TBoundaryConditions** boundaryConditions,TGeneEleMutationProbsLaw* geneEleMutationProbsLaw, TIndexGene mutationGene, TBoolean moreThanOnePointSubstitution);
		void DoGenePointSubstitution(TMutationLaw** mutationPDF,TBoundaryConditions** boundaryConditions,TGeneEleMutationProbsLaw* geneEleMutationProbsLaw, TBoolean moreThanOnePointSubstitution);
		void DoAllPointMutations(TSimulationParameters * simulationParameters);
		//void DoClustersDescriptionFromGenes();
		//static int orderFunctionPhenotype(const void* A_i, const void* B_i);
		void CopyOtherIndividual(CIndividual* other);
		void ComputeFitness(TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData ,  CFitnessEvaluator * fitnessComputer);//TFitnessMode fitnessMode, TGeneElement* objectiveFunction, TGenome genomeToEvaluate );
		TData* ComputePhenotype( CFitnessEvaluator * fitnessComputer);
		void ComputeCodingRatio();
		//void permute2PearlsElements(TGene perl1, TGene perl2,TIndexGeneElement cutInGene);
		//TIndexGene ComputeTransposonSize(TIndexGene cutGeneStart, TIndexGene cutGeneEnd, TIndexGene len);
		//TIndexGene PositionAfterTransformation(TIndexGene geneInTransposon, TIndexGene insertionGene, TIndexGene cutGeneStart,TIndexGene len);
		//void RearrangementMutationInduced(TIndexGene len, TIndexGene transposonSize, TIndexGene cutGeneStart,TIndexGene  insertionGene, TIndexGene cutGeneEnd);
		void MakeTranslocationPermutation(TIndexGene cutGeneStart, TIndexGene cutGeneEnd,  TIndexGene insertionGene, TIndexGene len, TIndexGene transposonSize,TIntergenicCut* intergenicCutParameter);
		TIndexGene ComputeTransposonSize(TIndexGene cutGeneStart, TIndexGene cutGeneEnd, TIndexGene len);
		void InjectivePermutationGeneElements(TGene* pearl1, TGene* pearl2, TIndexGeneElement cut, TGene* abstractGenesOrder);
		void BijectivePermutationGeneElements(TGene* pearl1, TGene* pearl2, TIndexGeneElement cut, TGene* abstractGenesOrder);
		void MakeDeletionPermutation(TIndexGene cutGeneStart, TIndexGene cutGeneEnd, TIndexGene len, TIndexGene transposonSize,TIntergenicCut* intergenicCutParameter);
		void MakeDuplicationPermutation(TIndexGene cutGeneStart, TGene* cutGeneEndInTransposon, TIndexGene insertionGene,  TIndexGene len, TIndexGene transposonSize,TIntergenicCut* intergenicCutParameter);
		void ChangeGenomeLength(TIndexGene newlength);
		void SetGeneElement(TIndexGene geneIndex, TIndexGeneElement elementIndex, TGeneElement value);
		void Shuffle_genome();
		void SaveGenome(TFileName nameSaveFile);
		void SaveFeatures(TFileName nameSaveFile,TGenerationsIndex curentGenerationIndex);
		void SavePhenotype(CFitnessEvaluator * fitnessComputer, TFileName nameSaveFile);
		void SaveDataSet(TData* dataset, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData, TFileName nameSaveFile);
		void SaveClassifiedData(TData* dataStream, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData ,CFitnessEvaluator * fitnessComputer, TFileName nameSaveFile);
		void ModifyPseudogenes(TMutationLaw** newPseudogenePDF,TBoundaryConditions** geneBoundaryConditions);
		static int  OrderFunctionGenotype(const void* A_i, const void* B_i);
		void OrderGenome();
};
 
#endif


