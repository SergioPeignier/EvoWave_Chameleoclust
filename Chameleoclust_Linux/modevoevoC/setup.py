#!/usr/bin/env python
import os
os.environ["CC"] = "c++"
#os.environ["CC"] = "nvcc"
from distutils.core import setup, Extension
module = Extension('modevoevo_c', ["modevoevo_c.cpp", "CPrng.cpp", "CIndividual.cpp", "CPopulation.cpp", "CSimulation.cpp",  "CFitnessEvaluator.cpp", "CFitnessModeTest1.cpp", "CFitnessCluster.cpp","CFitnessClusterKmeans.cpp","CFitnessRandom.cpp"],libraries=['gsl','blas','z'])
module.extra_compile_args = [ '-std=gnu++0x']
setup(name='modevoevo_c',
	version='1.0',
	ext_modules=[module])

#c++ main.cpp CIndividual.cpp CPrng.cpp CFitnessEvaluator.cpp  CFitnessCluster.cpp -lgsl -lgslcblas -lm -std=gnu++0x -lz
