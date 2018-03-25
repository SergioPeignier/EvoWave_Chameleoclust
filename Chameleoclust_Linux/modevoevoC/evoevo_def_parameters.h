#include <iostream>
#include <string>
#include "evoevo_def_state_pop.h"
#ifndef PARAM
#define PARAM

#define FITNESS_MODE_TEST_1 0
#define FITNESS_MODE_TEST_2 1
#define FITNESS_MODE_TEST_3 2
#define FITNESS_MODE_CLUSTERING_KMEANS 3
#define FITNESS_MODE_CLUSTERING 4
#define FITNESS_RANDOM_MODE_TEST 5
#define LAW_NORMAL 0
#define LAW_UNIFORM 1
#define LAW_TRANSITION_MATRIX 2

typedef float         TMutationRate;
typedef float         Tseed;
typedef float         TNormExponent;
typedef float         TSelectionStrength;
typedef float         TDimensionWeight;
typedef float         TStatisticMomentum;
typedef float         TTransitionProbability;
typedef short         TIndexTransitionProbability; 
typedef int           TPopulationSize;
typedef unsigned int  TGenerationsIndex;
typedef unsigned char TFitnessMode;
typedef unsigned char TLawFamilly;
typedef unsigned char TKmeansIterations;
typedef char*         TLogFile;
typedef char*         TBackupFileName;


/**
* \struct TTransitionMatrix
* \brief Contains a transition matrix
*/
typedef struct TTransitionMatrix{
	TIndexTransitionProbability matrixSize;
	TTransitionProbability** matrix;
}TTransitionMatrix;

/**
* \struct TMutationLawParameters
* \brief Contains informations concerning a mutation law
*/
typedef struct TMutationLawParameters{
	TGeneElement min;
	TGeneElement max;
	TStatisticMomentum mean;
	TStatisticMomentum standardDeviation;
	TLawFamilly law; 
	TTransitionMatrix* p_transitionMatrix;
}TMutationLawParameters;

/**
* \struct TMutationLaw
* \brief Contains the TMutationLawParameters for each tuple element
*/
typedef struct  TMutationLaw{
	TMutationLawParameters * mutationLaw;
	}TMutationLaw;

/**
* \struct TBoundaryConditions
* \brief Contains a tuple element bounds 
*/
typedef struct TBoundaryConditions{
	TGeneElement min;
	TGeneElement max;
	}TBoundaryConditions;

/**
* \struct TIntergenicCut
* \brief Contains the parameters to do an intergenic cut operation
*/
typedef struct TIntergenicCut{
	TGene* genesAbstractOrder;
	TBoolean cut;
	TGeneElement mincutpoint;
	TGeneElement maxcutpoint;
	}TIntergenicCut;

/**
* \struct TGeneEleMutationProbsLaw
* \brief Contains two TMutationRate pointers, one containing the tuple elements mutation probabilities (probGeneElementMutation) and another to run a wheel of fortune choice 
*/
typedef struct TGeneEleMutationProbsLaw{
	TMutationRate* probGeneElementMutation;
	TMutationRate* probGeneElementMutationWheelOfFortune;
	}TGeneEleMutationProbsLaw;

/**
* \struct TSimulationParameters
* \brief Contains the parameters that are required by the Simulation class
*/
typedef struct TSimulationParameters{
	Tseed                     prngSeed;
	TLogFile                  logFile;
	
	TIndexGeneElement         geneSize;
	TIndexGene                numberOfInitialGenes;
	TIndexPoint               sizeDataBuffer;
	TIndexCoordinate          sizeDataPoint;
	TIndexGene                genomeSizeLimit;
	
	TMutationLaw**            mutationLaws;                
    TMutationLaw**            initPositionsPDF;
    TBoundaryConditions**     boundaryConditions;
    TGeneEleMutationProbsLaw* geneEleMutationProbsLaw;
    
    TMutationRate             largeDuplicationRate;
    TMutationRate             largeDeletionRate;
    TMutationRate             largeTranslocationRate;
    TMutationRate             largeInvertionRate;
    TMutationRate             pointSubstitutionRate;      
    TMutationRate             pointInsertionRate;        
    TMutationRate             pointDeletionRate;          
    TBoolean                  moreThanOnePointSubstitution; 

	TPopulationSize           populationSize;
	TGenerationsIndex         currentGenerationIndex;
	TSelectionStrength        selectionPressure;

	TKmeansIterations         kmeansIterations;
	TNormExponent             normExponent;
	TBoolean                  sortByPosition;
	TIntergenicCut*           intergenicCutParameter;
	TBoolean                  saveBestIndividual;
	TFitnessMode              fitnessMode;
	TBoolean                  shuffleGenome;
    TFitness                  nonCodingGenomeFitness;
    TMutationLaw*             unknownDimentionRandomValuesGenerator;
    
    TDimensionWeight          residualDistanceWeight;
    
} TSimulationParameters;


#endif


