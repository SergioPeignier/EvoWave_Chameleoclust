#!/usr/bin/python
import sys
import time
import os
import pickle
import numpy as np
from definitions import *

def write_txt_to_file(filename_i,txt,rewrite = False):
	if rewrite: self.erase_file_content(filename_i)
	file_to_save = open(filename_i, "a")
	if isinstance(txt,str):
		file_to_save.write(txt)
	if isinstance(txt,list) or isinstance(txt,tuple):
		txt = np.array(txt)
	if isinstance(txt, (np.ndarray, np.generic) ):
		np.savetxt(file_to_save, txt, delimiter=SEP_IN_FILE, newline=END_OF_LINE_FILE, header='', footer='', comments=COMMENTS_IN_FILE)
	file_to_save.close()

def erase_file_content(file_name):
	file_to_save = open(file_name,"w")
	file_to_save.write("")
	file_to_save.close()

def save_object_to_file(object,filename):
	f = open(filename, 'wb')
	pickle.dump(object,f)
	f.close()
	
def load_object_from_file(filename):
	f      = open(filename,"rb")
	object = pickle.load(f)
	f.close()
	return object	

def load_simple_float_matrix_from_txt(filename):
	f = open(filename, "r")
	ans = [[float(field) for field in line.split(",")[0:-1]] for line in f]
	f.close()
	return ans	