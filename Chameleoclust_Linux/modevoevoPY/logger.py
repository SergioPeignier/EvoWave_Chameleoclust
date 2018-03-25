#!/usr/bin/python
import time
import os
import resource
from definitions import *
from saveloadfunctions import * 

class logger:
	def __init__(self,exp_params_title={},parent_folder_name_params={},parent_folder_name = "./",folder_name="./",time_name = "time_generations_evoevo.txt"):
		self.parent_folder_name = parent_folder_name
		self.parent_folder_name_params = parent_folder_name_params
		self.folder_name = folder_name
		self.time_name = time_name
		self.time0 = time.time()
		self.exp_params_title = exp_params_title
		self.exp_params_title["TIME"] = time.strftime("%d_%m_%Y_%X")
		if self.folder_name == "": self.folder_name = self.generatefoldername(self.exp_params_title, BASAL_NAME_FOLDER)
		if self.parent_folder_name == "": self.parent_folder_name = self.generatefoldername(self.parent_folder_name_params,BASAL_NAME_PARENT_FOLDER)
		self.path_to_dir(os.path.join(self.parent_folder_name,self.folder_name))
	
	def generatefoldername(self,dico,basalname = BASAL_NAME_FOLDER):
		return "_".join([basalname]+["_".join(map(str,x)) for x in dico.items()])
	
	def path_to_dir(self,path):
		if not os.path.exists(path):
			l=[]
			p = "/"
			l = path.split("/")
			i = 1
			p = l[0]+"/"
			while i < len(l):
				p = p + l[i] + "/"
				i = i + 1
				if not os.path.exists(p):
					os.mkdir(p, 0755)

	def save_object_to_file(self,object,filename):
		save_object_to_file(object,self.get_full_name(filename))
	
	def load_object_from_file(self,filename):
		return load_object_from_file(self.get_full_name(filename))
		
	def get_full_name(self,filename_i):
		return os.path.join(self.parent_folder_name,self.folder_name,filename_i)
		
	def write_time_mesures(self,generation):
		usage=resource.getrusage(resource.RUSAGE_SELF)
		txt_line = [[generation,usage[0],usage[1],(usage[2]*resource.getpagesize())/1000000.0]]
		write_txt_to_file(self.get_full_name(self.time_name),txt_line,False)
"""
def generatestandardlogger(datasetname, eparams,paramsexp = {},paramsparent={}):
	exp_params_title={"seed": eparams.prng_seed}
	exp_params_title.update(paramsexp)
	parent_folder_name_params = {"dataset":datasetname}
	parent_folder_name_params.update(paramsparent)
	return logger(exp_params_title=exp_params_title, parent_folder_name_params= parent_folder_name_params)
"""
