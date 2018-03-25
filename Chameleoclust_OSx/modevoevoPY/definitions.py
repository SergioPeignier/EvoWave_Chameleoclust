BASAL_NAME_FOLDER   = "RUN_EVOEVO"
BASAL_NAME_PARENT_FOLDER   = "./RESULTS"
FILE_NAME_SEPARATOR = "_"
SEP_IN_FILE         = "\t"
END_OF_LINE_FILE    = "\n"
COMMENTS_IN_FILE    = '# '

GAUSSIAN            =   0
UNIFORM             =   1
TRANSITION_MATRIX   =   2
GENE_SIZE           =   4
MAX_GENE_SIZE       =   4

FITNESS_COPY_FUNCTION_WITH_DOTS   = 0
FITNESS_COPY_FUNCTION_TRIANGLES   = 1
FITNESS_AEVOL                     = 2
FITNESS_CLUSTER                   = 4
FITNESS_RANDOM                    = 5

PSEUDOGENE_GENE_TYPE = 0
CENTROID_GENE_TYPE   = 1
DIMENSION_GENE_TYPE  = 2
POSITION_GENE_TYPE   = 3


MAX_INT = 100000000
MIN_INT = -100000000

STATS_FITNESS           = 0
STATS_CODING_RATIO      = 1
STATS_GENOME_LENGTH     = 2
STATS_FITNESS_INDEX_POS = 3
STATS_PARENT            = 4
STATS_INDEX             = 5

STATS_FEATURES = ["fitness",
				  "coding_ratio",
				  "genome_length",
				  "fitness_pos",
				  "parent_index",
				  "index"]
				  
STATS_TO_DROP  = ["fitness_pos"]

STD_AGGREGATION_COLUMNS =  ["fitness",
							"coding_ratio",
							"genome_length",
							"best_fitness",
							"best_coding_ratio",
							"best_genome_length",
							"best_nb_core_points",
							"best_avg_core_points_dim",
							"best_nb_clusters",
							"best_avg_clusters_dim",
							"entropy",
							"accuracy",
							"f1",
							"CE"]
							
STD_USEFULL_STATS_AGGREGATION = ["fitness",
								 "coding_ratio",
								 "genome_length"]
								 
STD_FUNCTIONAL_EVALUATION_COLUMNS = ["Entropy",
									 "Accuracy",
									 "F1",
									 "CE",
									 "RNIA",
									 "CE_SSC",
									 "Coverage"]
									 
STD_STRUCTURAL_EVALUATION_COLUMNS = ["nb_core_points",
									 "core_points_avg_dim",
									 "nb_clusters",
									 "clusters_avg_dim"]
									 
STD_FITNESS_RANDOM_WALK_STUDY = ["rho_f","tau_f"]

STD_Fw_RANDOM_WALK_STUDY      = ["rho_fw","tau_fw"]

STD_Fv_RANDOM_WALK_STUDY      = ["rho_fv","tau_fv"]

STD_Fb_RANDOM_WALK_STUDY      = ["rho_fb","tau_fb"]

STD_FW_FV_FB_COLUMNS          = ["fw","fv","fb"]

NONABSORBINGSTATE   = [[0,1],
                       [1,0]]
                       
PSEUDOGENEABSORBING = [[1,0],
                       [1,0]]
                       
GENEABSORBING       = [[0,1],
                       [0,1]]
                       
BOTHABSORBING       = [[1,0],
                       [0,1]]
                       
STD_CURRENT_GENERATION_INDEX  = 0
STD_PRNG_SEED                 = 0                   
STD_MUT_RATE                  = 0.00142
STD_GENE_SIZE                 = 4
STD_INIT_GENES                = 200
STD_DATA_BUFFER               = 1024
STD_SIZE_DATA_POINT           = 1024
STD_GENOME_LIM_SIZE           = 50000
STD_POP_SIZE                  = 300
STD_SELECT_PRESSURE           = 0.5
STD_KMEAN_ITERATIONS          = 0
STD_NORM_EXPONENT             = 1
STD_FITNESS_MODE              = FITNESS_CLUSTER 
STD_LOG_FILE                  = "test.PICKLE"


STD_NON_CODING_GENOME_FITNESS = -100000000
STD_SHUFFLE_GENOME            = False
STD_CASCADE_MUT               = False
STD_SAVE_BEST_INDIVIDUAL      = True
STD_SORT_BY_POSITION          = True

STD_WINDOW_SIZE               = 0.1

STD_GENE_MUTATION_PROBABILITIES = [0.25,0.25,0.25,0.25]

STD_ORDER_FOR_INTERGENIC_CUT_OFF = [0,1,2,3]

STD_PRECISION = 1000



