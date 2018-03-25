#ifndef STATE
#define STATE

#define GENE_SIZE 5
#define INDEX_TYPE 0 
#define INDEX_CLUST 1 
#define INDEX_DIM 2 
#define INDEX_POS 3 
#define INDEX_WEIG 4
#define CENTROID_GENE_TYPE 1
#define PSEUDOGENE_GENE_TYPE 0
#define SIZE_TEMP_GENE_ARRAY 5000
#define COORDINATE_LENGTH 2 
#define DIMENSION_IN_COORDINATE 0
#define POSITION_IN_COORDINATE 1
#define NB_INDIVIDUAL_FEATURES 6
#define FEATURE_FITNESS_POS 0
#define FEATURE_CODING_RATIO_POS 1
#define FEATURE_GENOME_LENGTH_POS 2
#define FEATURE_FITNESS_INDEX_POS 3
#define FEATURE_PARENT_INDEX_POS 4
#define FEATURE_INDIV_INDEX_POS 5
#define NB_LARGE_MUTATIONS 4
#define LARGE_DUPLICATIONS_NB 0 
#define LARGE_DELETIONS_NB 1 
#define LARGE_TRANSLOCATIONS_NB 2
#define LARGE_INVERTIONS_NB 3 

typedef int TCoordinateElement;
typedef int TIndexCoordinateElement;
typedef int TIndexPoint; 
typedef int TIndexCoordinate;
typedef int TIndexIndividual;
typedef int TIndexGene;
typedef int TIndexSimulations;
typedef int TIndexGeneElement;
typedef double TReproductionProb;
typedef float TFitness;
typedef float TCodingRatio;
typedef int TGeneElement; 

class CIndividual;
class CPopulation;

/**
* \struct TSimulationState
* \brief Contains the simulation population (CPopulation *)
*/
typedef struct TSimulationState{
	CPopulation * r_population;
} TSimulationState;

/**
* \struct TArrayIndividuals
* \brief  Contains an array of pointers to CIndividual class instances and their respective indexes
*/
typedef struct TArrayIndividuals{
	TIndexIndividual nbElts; // Number of elements actually used in the array
	CIndividual ** p_elt; // array[MAX_NB_INDIVIDUAL]
} TArrayIndividuals;

/**
* \struct TGene
* \brief  TGeneElement array of size GENE_SIZE : a genome tuple
*/
typedef struct TGene{
	TGeneElement gene[GENE_SIZE];
	}TGene;

/**
* \struct TArrayGenes
* \brief  Array of TGenes, number of elements and maximal number of elements that it can contain 
*/
typedef struct TArrayGenes{
    TIndexGene nbElts;
    TIndexGene maxNbElts;
    TGene * p_elt;				
} TArrayGenes;

/**
* \struct TGenome
* \brief  Pointer to TArrayGenes
*/
typedef struct TGenome{
    TArrayGenes * p_arrayGenes;
}TGenome;

/**
* \struct TCoordinate
* \brief  Array of coordinate elements of size 2
*/
typedef struct TCoordinate{
	TCoordinateElement coordinate[COORDINATE_LENGTH];
}TCoordinate;

/**
* \struct TPointCoordinatesArray
* \brief  Array of TCoordinate structs and a the number of element it contains
*/
typedef struct TPointCoordinatesArray{
    TIndexCoordinate nbElts;
    TCoordinate * p_elt;	
}TPointCoordinatesArray;

/**
* \struct TPoint
* \brief  Array of TPointCoordinatesArray structs, a label and a cluster id
*/
typedef struct TPoint{
	TPointCoordinatesArray * p_coordinates;
	TGeneElement label;
	TGeneElement cluster;
}TPoint;

/**
* \struct TArrayPoints
* \brief  Array of TPoint structs and the number of elements it contains
*/
typedef struct TArrayPoints{
    TIndexPoint nbElts;
    TPoint * p_elt;				
} TArrayPoints;

/**
* \struct TData
* \brief  pointer to a TArrayPoints struct
*/
typedef struct TData{
    TArrayPoints * p_arrayPoints;
}TData;

/**
* \struct TDimAggregator
* \brief  TCoordinateElement array, the number of elements it contains and the dimension it deals with
*/
typedef struct TDimAggregator{
    TIndexCoordinate nbElts;
    TCoordinateElement * p_elt;	
    TCoordinateElement dimension;
}TDimAggregator;

/**
* \struct TPointAggregator
* \brief  Array of TDimAggregators, the number of elements it contains and the id
*/	
typedef struct TPointAggregator{
	TDimAggregator * dimsArray;
	TIndexCoordinate nbElts;
	TGeneElement id;
}TPointAggregator;

/**
* \struct TArrayPointsAggregator
* \brief  Array of TPointAggregator, the number of elements it contains
*/	
typedef struct TArrayPointsAggregator{
    TIndexPoint nbElts;
    TPointAggregator * p_elt;				
} TArrayPointsAggregator;

/**
* \struct TDataAggregator
* \brief Pointer of TArrayPointsAggregator
*/
typedef struct TDataAggregator{
    TArrayPointsAggregator * p_arrayPoints;
}TDataAggregator;
	
/**
* \struct TIndividualFeatures
* \brief Individual fitness, functional ratio, genome size, rank, parent index, index
*/
typedef struct TIndividualFeatures{
	TFitness           fitness;
	TCodingRatio       codingRatio;
	TIndexGene         genomeSize;
	TIndexIndividual   fitnessIndex; //may be erased latter
	TIndexIndividual   parentIndex;
	TIndexIndividual   individualIndex;
}TIndividualFeatures;

#endif
