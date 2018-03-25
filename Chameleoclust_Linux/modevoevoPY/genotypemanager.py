#!/usr/bin/python
from egenomeparameters import *
from saveloadfunctions import *


def eloadgenomefromfile(parentalgenomefile,extension="txt"):
	genome = []
	if extension == "txt":
		genome = load_simple_float_matrix_from_txt(parentalgenomefile)
	if extension == "PICKLE":
		genome = load_object_from_file(parentalgenomefile)
	return genome	

def egeneraterandomDNAstrand(length,pseudognesgenerationlaw):
	return([[law.drawrandomnumber() for law in pseudognesgenerationlaw] for i in range(length)])

def econverttopseudogenes(genome):
	for i,gen in enumerate(genome):
		genome[i][PSEUDOGENE_GENE_TYPE] = 0
	return genome

def egeneraterandomDNAfromindividualDNAdistribution(length,genome):
	composition = zip(*genome)
	return [ [random.sample(c,1)[0] for c in composition] for i in range(length) ]

def ecomparegenomes(genome1,genome2,only_coding=0):
	N = 0
	for gene1 in genome1:
		for gene2 in genome2:
			if not only_coding:
				i = 0
				equalgenes = 1
				while i < len(gene1):
					if gene1[i] != gene2[i]:
						equalgenes = 0
						break
					i += 1
				N += equalgenes
			if only_coding and gene1[PSEUDOGENE_GENE_TYPE] and gene2[PSEUDOGENE_GENE_TYPE]:
				i = 0
				equalgenes = 1
				while i < len(gene1):
					if gene1[i] != gene2[i]:
						equalgenes = 0
						break
					i += 1
				N += equalgenes
	return N

def ecomputecodinglength(genome):
	c = 0
	for gene in genome:
		if gene[PSEUDOGENE_GENE_TYPE]:
			c += 1
	return c


def eerasetuplewithtypex(genome,type_tuples_to_erase = 0,nb_tuples_to_erase="all"):
	new_genome	 = []
	tuples_with_type = []
	for i,gen in enumerate(genome):
		if gen[PSEUDOGENE_GENE_TYPE] == type_tuples_to_erase:
			tuples_with_type.append(i)
	if len(tuples_with_type) <= nb_tuples_to_erase or nb_tuples_to_erase=="all":
		nb_tuples_to_erase = len(tuples_with_type)
	to_erase = random.sample(tuples_with_type, nb_tuples_to_erase)
	for i,gene in enumerate(genome):
		if i not in to_erase:
			new_genome.append(gene[:])
	return new_genome

def ecompute_nb_genes(genome):
	nb_genes = 0
	for gene in genome:
		if gene[PSEUDOGENE_GENE_TYPE]:
			nb_genes += 1
	return nb_genes

def ecompute_nb_pseudogenes(genome):
	nb_genes = 0
	for gene in genome:
		if not gene[PSEUDOGENE_GENE_TYPE]:
			nb_genes += 1
	return nb_genes

def ereducesize(genome,toreduce):
	genome = random.sample(genome,len(genome))
	genome = genome[toreduce:len(genome)]
	return genome

def eerasepseudogenes(genome,nb_pseudogenes_to_erase="all"):
	return eerasetuplewithtypex(genome,nb_tuples_to_erase=nb_pseudogenes_to_erase)

def eerasegenes(genome,nb_genes_to_erase="all"):
	return eerasetuplewithtypex(genome,type_tuples_to_erase = 1,nb_tuples_to_erase=nb_genes_to_erase)

def egeneraterandomDNAfromindividualNonCodingDNA(length,individualgenome):
	noncodingcontent = eerasegenes(individualgenome)
	composition = zip(*noncodingcontent)
	return [ [random.sample(c,1)[0] for c in composition] for i in range(length) ]

