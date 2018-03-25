#include "evoevo_def_state_pop.h"
#include "evoevo_def_base.h"
#include "evoevo_def_parameters.h"
#include "CFitnessEvaluator.h"
#include "CPrng.h"

#ifndef FITNESSCLUSTERKmeans
#define FITNESSCLUSTERKmeans
class CFitnessClusterKmeans : public CFitnessEvaluator{
	public:
		TGenome*       a_codingGenes;
		TMutationLaw*  a_statisticalLaw;
		TNormExponent  a_norm;
		TFitness       a_nonCodingGenomeFitness;
		TIndexGene     a_start;
		TIndexGene     a_stop;
		TDataAggregator*   a_modified_phenotype;
		TPoint*        a_phenotypeModifier;
		TPoint*        a_phenotypeModifierlocal;
		TBoolean       a_convergence_k_means;
		TKmeansIterations a_KMeanIterations;
		TGeneElement   a_maxNbDims;
		TIndexPoint    a_maxNbPoints;
		CFitnessClusterKmeans(CPrng* prng,TIndexGene maxNbGenes, TFitness nonCodingGenomeFitness, TNormExponent norm, TMutationLaw* statisticalLaw, TKmeansIterations nbMaxKmeansIterations,TGeneElement maxNbClust,  TGeneElement maxNbDims,  TGeneElement maxNbPoints);
		~CFitnessClusterKmeans ();
		TFitness EvaluateFitness(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		void AsociatePleiotropicGenes(TGenome* genomeCoding);
		TGeneElement AsociateGenesWeightedSum(TGenome* genomeCoding, TIndexGene start,TIndexGene  stop,TIndexGeneElement geneSize);
		void CutAndSavePartOfGenome(TGenome* genomeSorted, TIndexGene start, TIndexGene stop);
		void GetBoundaryGenesOfTypeX(TGenome* genomeSorted, TGeneElement type);
		bool HaveSameType_Cluster_Dimention(TGene gene1, TGene gene2);
		static int orderFunctionPhenotype(const void* A_i, const void* B_i);
		void TranscriptionGenome2CodingGenome(TGenome* genome, TGenome* genomeCoding);
		void PrintGenome(TGenome* genome);
		TFitness AddtoDistance(TFitness additional, TNormExponent norm);
		TFitness AddtoDistanceNotInData();
		TFitness DistanceNormalization( TFitness distance ,TGeneElement dimensionality,TNormExponent  norm);
		TData* ComputePhenotype(TGenome* genomeToEvaluate,TIndexGeneElement geneSize);
		void ForgetPhenotype(TData* phenotype);
		TData* ClassifyData(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		TFitness EvaluateFitnessInDataSubset(TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		void PrintPoint(TPoint point);
		void Coding2Clusters(TGenome* genomeCoding,TData* clusters,TIndexGeneElement geneSize);
		TFitness ComputeAllDistancesTo0(TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		void PrintDataSet(TData* a_dataStream,TIndexPoint a_startObservationPosInData,TIndexPoint a_stopObservationPosInData);
		void UpdatePhenotype();
		TFitness AssignPoints(TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		void InitializePhenotypeModified();
		void PrintDataSetToFile(FILE* f,TData* a_dataStream,TIndexPoint a_startObservationPosInData,TIndexPoint a_stopObservationPosInData);
		void AllocateDataAggregatorMem(TDataAggregator* dataAggregator, TGeneElement maxNbClust,  TGeneElement maxNbDims,  TGeneElement maxNbPoints);
		void freeDataAggregator(TDataAggregator* dataAggregator);
		TCoordinateElement qselect(TCoordinateElement *arr, TIndexPoint n, TIndexPoint k);
};
#endif
