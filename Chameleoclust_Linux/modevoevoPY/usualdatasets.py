from dataset import *
from aevol_data_generator import *
breast = {"datatable_file_name" : "../../../data/Databases/real_world_data/breast.arff.txt",
"hidden_cluster_index" : [33],
"true_file_name" : "../../../data/Databases/real_world_data/breast.true"}

diabetes = {"datatable_file_name" : "../../../data/Databases/real_world_data/diabetes.arff.txt",
"hidden_cluster_index" : [8],
"true_file_name" : "../../../data/Databases/real_world_data/diabetes.true"}

glass = {"datatable_file_name" : "../../../data/Databases/real_world_data/glass.arff.txt",
"hidden_cluster_index" : [9],
"true_file_name" : "../../../data/Databases/real_world_data/glass.true"}

liver = {"datatable_file_name" : "../../../data/Databases/real_world_data/liver.arff.txt",
"hidden_cluster_index" : [6],
"true_file_name" : "../../../data/Databases/real_world_data/liver.true"}

pendigits = {"datatable_file_name" : "../../../data/Databases/real_world_data/pendigits.arff.txt",
"hidden_cluster_index" : [16],
"true_file_name" : "../../../data/Databases/real_world_data/pendigits.true"}

shape = {"datatable_file_name" : "../../../data/Databases/real_world_data/shape.arff.txt",
"hidden_cluster_index" : [17],
"true_file_name" : "../../../data/Databases/real_world_data/shape.true"}

vowel = {"datatable_file_name" : "../../../data/Databases/real_world_data/vowel.arff.txt",
"hidden_cluster_index" : [10],
"true_file_name" : "../../../data/Databases/real_world_data/vowel.true"}

S1500 = {"datatable_file_name" : "../../../data/Databases/synth_dbsizescale/S1500.arff",
"hidden_cluster_index" : [20],
"true_file_name" : "../../../data/Databases/synth_dbsizescale/S1500.true"}

S2500 = {"datatable_file_name" : "../../../data/Databases/synth_dbsizescale/S2500.arff",
"hidden_cluster_index" : [20],
"true_file_name" : "../../../data/Databases/synth_dbsizescale/S2500.true"} 

S3500 =  {"datatable_file_name" : "../../../data/Databases/synth_dbsizescale/S3500.arff",
"hidden_cluster_index" : [20],
"true_file_name" : "../../../data/Databases/synth_dbsizescale/S3500.true"}

S4500 =  {"datatable_file_name" : "../../../data/Databases/synth_dbsizescale/S4500.arff",
"hidden_cluster_index" : [20],
"true_file_name" : "../../../data/Databases/synth_dbsizescale/S4500.true"}

S5500 =  {"datatable_file_name" : "../../../data/Databases/synth_dbsizescale/S5500.arff",
"hidden_cluster_index" : [20],
"true_file_name" : "../../../data/Databases/synth_dbsizescale/S5500.true"}

D05 =  {"datatable_file_name" : "../../../data/Databases/synth_dimscale/D05.arff",
"hidden_cluster_index" : [5],
"true_file_name" : "../../../data/Databases/synth_dimscale/D05.true"}

D10 =  {"datatable_file_name" : "../../../data/Databases/synth_dimscale/D10.arff",
"hidden_cluster_index" : [10],
"true_file_name" : "../../../data/Databases/synth_dimscale/D10.true"}

D15 =  {"datatable_file_name" : "../../../data/Databases/synth_dimscale/D15.arff",
"hidden_cluster_index" : [15],
"true_file_name" : "../../../data/Databases/synth_dimscale/D15.true"}

D20 =  {"datatable_file_name" : "../../../data/Databases/synth_dimscale/D20.arff",
"hidden_cluster_index" : [20],
"true_file_name" : "../../../data/Databases/synth_dimscale/D20.true"}

D25 =  {"datatable_file_name" : "../../../data/Databases/synth_dimscale/D25.arff",
"hidden_cluster_index" : [25],
"true_file_name" : "../../../data/Databases/synth_dimscale/D25.true"}

D50 =  {"datatable_file_name" : "../../../data/Databases/synth_dimscale/D50.arff",
"hidden_cluster_index" : [50],
"true_file_name" : "../../../data/Databases/synth_dimscale/D50.true"}

D75 =  {"datatable_file_name" : "../../../data/Databases/synth_dimscale/D75.arff",
"hidden_cluster_index" : [75],
"true_file_name" : "../../../data/Databases/synth_dimscale/D75.true"}

N10 =  {"datatable_file_name" : "../../../data/Databases/synth_noisescale/N10.arff",
"hidden_cluster_index" : [20],
"true_file_name" : "../../../data/Databases/synth_noisescale/N10.true"}

N30 =  {"datatable_file_name" : "../../../data/Databases/synth_noisescale/N30.arff",
"hidden_cluster_index" : [20],
"true_file_name" : "../../../data/Databases/synth_noisescale/N30.true"}

N50 =  {"datatable_file_name" : "../../../data/Databases/synth_noisescale/N50.arff",
"hidden_cluster_index" : [20],
"true_file_name" : "../../../data/Databases/synth_noisescale/N50.true"}

N70 =  {"datatable_file_name" : "../../../data/Databases/synth_noisescale/N70.arff",
"hidden_cluster_index" : [20],
"true_file_name" : "../../../data/Databases/synth_noisescale/N70.true"}


real_synth_data_sets = {"breast" : breast,
						"diabetes": diabetes,
						"glass": glass,
						"liver": liver,
						"pendigits": pendigits,
						"shape": shape,
						"vowel": vowel,
						"S1500": S1500,
						"S2500": S2500,
						"S3500": S3500,
						"S4500": S4500,
						"S5500": S5500,
						"D05": D05,
						"D10": D10,
						"D15": D15,
						"D20": D20,
						"D25": D25,
						"D50": D50,
						"D75": D75,
						"N10": N10,
						"N30": N30,
						"N50": N50,
						"N70": N70}

def get_usual_datasets(datatable_file_name,hidden_cluster_index,found_cluster = [],comments = "@",extension = "arff", sep = ",", true_file_name = ""):
	a = data(datatable_file_name,
			 hidden_cluster=hidden_cluster_index,
			 found_cluster=found_cluster,
			 comments = comments,
			 extension = extension,
			 sep = sep)
	a.standardize_table()
	a.multiply_table_by_precision_and_convert_to_int()
	a.table_to_data()
	tf = []
	if true_file_name != "":
		tf = true_file(true_file_name,len(a.features))
		tf.generate_dico_arff_true_clusters(a.hidden_cluster_membership)
		tf.compute_df_true_dico(a.conversion_dico[a.hidden_cluster[0]])
	return {"arff": a, "true": tf}
	
	
get_breast = lambda : get_usual_datasets(** real_synth_data_sets["breast"])
get_diabetes = lambda : get_usual_datasets(** real_synth_data_sets["diabetes"])
get_glass = lambda : get_usual_datasets(** real_synth_data_sets["glass"])
get_liver = lambda : get_usual_datasets(** real_synth_data_sets["liver"])
get_pendigits = lambda : get_usual_datasets(** real_synth_data_sets["pendigits"])
get_shape = lambda : get_usual_datasets(** real_synth_data_sets["shape"])
get_vowel = lambda : get_usual_datasets(** real_synth_data_sets["vowel"])
get_pendigits = lambda : get_usual_datasets(** real_synth_data_sets["pendigits"])
get_shape = lambda : get_usual_datasets(** real_synth_data_sets["shape"])
get_vowel = lambda : get_usual_datasets(** real_synth_data_sets["vowel"])
get_S1500 = lambda : get_usual_datasets(** real_synth_data_sets["S1500"])
get_S2500 = lambda : get_usual_datasets(** real_synth_data_sets["S2500"])
get_S3500 = lambda : get_usual_datasets(** real_synth_data_sets["S3500"])
get_S4500 = lambda : get_usual_datasets(** real_synth_data_sets["S4500"])
get_S5500 = lambda : get_usual_datasets(** real_synth_data_sets["S5500"])
get_D05 = lambda : get_usual_datasets(** real_synth_data_sets["D05"])
get_D10 = lambda : get_usual_datasets(** real_synth_data_sets["D10"])
get_D15 = lambda : get_usual_datasets(** real_synth_data_sets["D15"])
get_D20 = lambda : get_usual_datasets(** real_synth_data_sets["D20"])
get_D25 = lambda : get_usual_datasets(** real_synth_data_sets["D25"])
get_D50 = lambda : get_usual_datasets(** real_synth_data_sets["D50"])
get_D75 = lambda : get_usual_datasets(** real_synth_data_sets["D75"])
get_N10 = lambda : get_usual_datasets(** real_synth_data_sets["N10"])
get_N30 = lambda : get_usual_datasets(** real_synth_data_sets["N30"])
get_N50 = lambda : get_usual_datasets(** real_synth_data_sets["N50"])
get_N70 = lambda : get_usual_datasets(** real_synth_data_sets["N70"])

get_aevol = lambda : {"arff": aevol_data(), "true": []}


def get_birds():
	a = data("../../../data/Databases/birds/pc_data.bird.txt",
		hidden_cluster=["0"],
		extension = "csv",
		sep = " ",
		index=None,
		header= 0)
	a.standardize_table()
	a.multiply_table_by_precision_and_convert_to_int()
	a.table_to_data()
	return {"arff": a, "true":[]}
	
def get_chemical():
	a = data("../../../data/Databases/ChameleoChem/descriptors.csv",
		hidden_cluster=[],
		extension = "csv",
		sep = ",",
		index=None,
		header= 0)
	a.multiply_table_by_precision_and_convert_to_int()
	a.table_to_data()
	return {"arff": a, "true":[]}

	
def get_bacterial(generation=""):
	if isinstance(generation,str):
		name = "../../../data/Databases/Bactos/clustering_data_total.txt"
	else:
		name = "../../../data/Databases/Bactos/clustering_data_"+str(generation)+".txt"
	a = data(name,
		useless = ["id","trophic_level","t"],
		hidden_cluster=["trophic_group"],
		extension = "csv",
		sep = " ",
		index=None,
		header= 0)
	a.standardize_table()
	a.multiply_table_by_precision_and_convert_to_int()
	a.table_to_data()
	return {"arff": a, "true":[]}



def get_birds_new(sex="male"):
	if sex == "female":
		name = "../../../data/Databases/birds/com.female.distance.calls.data"
	if sex == "male":
		name = "../../../data/Databases/birds/com.male.distance.calls.data"
	a = data(name,
		useless = range(21,100),
		hidden_cluster=[0],
		extension = "csv",
		sep = " ",
		index=None,
		header= None)
	a.standardize_table()
	a.multiply_table_by_precision_and_convert_to_int()
	a.table_to_data()
	return {"arff": a, "true":[]}
	
def get_chemical_T():
	a = data("../../../data/Databases/ChameleoChem/descriptors_T.csv",
		hidden_cluster=[],
		extension = "csv",
		sep = ",",
		index=None,
		header= 0)
	a.multiply_table_by_precision_and_convert_to_int()
	a.table_to_data()
	return {"arff": a, "true":[]}
