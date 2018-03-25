#include "evoevo_def_state_pop.h"
#include "evoevo_def_base.h"
#include "evoevo_def_parameters.h"
#include "CFitnessEvaluator.h"
#include "CPrng.h"

#ifndef FITNESSCLUSTER
#define FITNESSCLUSTER
class CFitnessCluster : public CFitnessEvaluator{
	public:
		TGenome*       a_codingGenes;
		TMutationLaw*  a_statisticalLaw;
		TNormExponent  a_norm;
		TFitness       a_nonCodingGenomeFitness;
		TIndexGene     a_start;
		TIndexGene     a_stop;
		CFitnessCluster(CPrng* prng,TIndexGene maxNbGenes, TFitness nonCodingGenomeFitness, TNormExponent norm, TMutationLaw* statisticalLaw, TGeneElement maxNbClust,  TGeneElement maxNbDims);
		~CFitnessCluster ();
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
		void ForgetPhenotype();
		TData* ClassifyData(TGenome* genomeToEvaluate,TIndexGeneElement geneSize, TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		TFitness EvaluateFitnessInDataSubset(TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		void PrintPoint(TPoint point);
		void Coding2Clusters(TGenome* genomeCoding,TData* clusters,TIndexGeneElement geneSize);
		TFitness ComputeAllDistancesTo0(TData* objectiveFunction, TIndexPoint startObservationPosInData, TIndexPoint stopObservationPosInData);
		void PrintDataSet(TData* a_dataStream,TIndexPoint a_startObservationPosInData,TIndexPoint a_stopObservationPosInData);

};
#endif
