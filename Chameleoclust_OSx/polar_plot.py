#-*- coding: utf-8 -*-
# 12/02/2016 : Print polar plots of the given datas. Build for Evowave datas
MAXNBPOINTS = 100
MAXNBCLUSTERS = 10
import itertools

import numpy as np
import pandas as pd
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import seaborn as sns
import math
class PolarPlotter(object):
    """
    Usage :
    PolarPlotter(data, classes, ...).plot()

    2 possibilites (only_medians = True/False:
        - 1 plot by cluster (one line by capture)
        - 1 line by cluster (median point), (one color by cluster)
    """
    def __init__(self,max_=0,min_=0, only_medians=False, r_base=0.5, ncols=5, palette=list(itertools.permutations([50,150,200])) + list(itertools.permutations([0,100,200]))):
        self.win = pg.GraphicsWindow(title="Polar plot of evowave datas")
        self.win.resize(1000,600)
    
        self.class_col = []
        self.only_medians = only_medians

        self.min = min_ #Les points du petit cercle ont une valeur de self.min
        self.max = max_

        self.n_dim = 0 #nb de colonnes
        self.n_classes = 0

        self.r_base = r_base #rayon du petit cercle autour de la figure
        self.ncols  = ncols #nombre de colonnes d'images
        self.palette = palette

    def polar2cartesian(self, r_, dim):
    	r = np.append(r_,r_[0])
        theta = (2 * np.pi * dim) / self.n_dim
        x = (r + self.r_base) * np.cos(theta)
        y = (r + self.r_base) * np.sin(theta)
        return x,y

    def point2shape(self, vec):
        return self.polar2cartesian(np.array(vec)-self.min, np.append(np.arange(len(vec)),0))

    def draw_limit_circles(self, p):
        p.addItem(pg.PlotDataItem(x=[0,0],y=[0,0], pen=None, symbol='o', symbolBrush=None, symbolSize=[(self.min+self.r_base)*2, (self.max+self.r_base)*2], pxMode=False))


class PolarPlotterAllInOne(PolarPlotter):
	def __init__(self,max_ = 10, min_ = -10 ,maximal_nb_clusters = MAXNBCLUSTERS, maximal_nb_points_in_plot = MAXNBPOINTS, r_base=0.5, ncols=5, palette=[[int(255*c) for c in col] for col in sns.color_palette("hls", MAXNBCLUSTERS)]):
		PolarPlotter.__init__(self,max_,min_, only_medians=maximal_nb_points_in_plot == 0, r_base=r_base, ncols=ncols, palette=palette)
		self.plot = self.win.addPlot(title="Clusters")
		self.maximal_distance = self.max-self.min+self.r_base
		self.plot.setRange(xRange=[-self.maximal_distance,self.maximal_distance], yRange = [-self.maximal_distance,self.maximal_distance])
		self.data_points_plot = [pg.PlotDataItem(pen=pg.mkPen(color=(200, 200, 200,30))) for _ in xrange(maximal_nb_points_in_plot)]
		self.data_clusters_plot = [pg.PlotDataItem(pen=pg.mkPen(color=self.palette[i])) for i in xrange(maximal_nb_clusters)]
		self.data = np.array([])
		for point_plot in self.data_points_plot:
			 self.plot.addItem(point_plot)
		for cluster_plot in self.data_clusters_plot:
			 self.plot.addItem(cluster_plot)	 
		self.oldest_point_index = 0
		self.maximal_nb_points_in_plot = maximal_nb_points_in_plot
		self.maximal_nb_clusters = maximal_nb_clusters
		
	def plot_point(self,point):
		point_plot_modified = -1
		if self.maximal_nb_points_in_plot > 0:
			if self.data.shape[0] < self.maximal_nb_points_in_plot or len(self.data.shape)==1:
				if self.data.shape[0] == 0 : 
					self.data = np.array(point)
					point_plot_modified = 0
				else:
					self.data = np.vstack((self.data,np.array(point)))
					point_plot_modified = self.data.shape[0]-1
			else:
				self.data[self.oldest_point_index] = np.array(point)
				point_plot_modified = self.oldest_point_index
				self.oldest_point_index += 1
				self.oldest_point_index = self.oldest_point_index % self.maximal_nb_points_in_plot
			self.n_dim = max(self.n_dim,len(point))
		#self.max = self.data.max().max()
		#self.min = self.data.min().min()
		x,y = self.point2shape(point)
		self.data_points_plot[point_plot_modified].setData(x,y)
        	
	def plot_cluster(self,cluster_coordinates,cluster_id):
		x,y = self.point2shape(cluster_coordinates)
		self.data_clusters_plot[cluster_id].setData(x,y)
		
	def plot_clusters(self,clusters_coordinates):
		for cluster_plot in self.data_clusters_plot:
			cluster_plot.setData([],[])
		for key in clusters_coordinates:
			self.plot_cluster(clusters_coordinates[key],key)
	
	def plot_dataset_points(self,dataset_points):
		for point in dataset_points:
			self.plot_point(point)

class PolarPlotterAllClustersInDifferentPlots(PolarPlotter):
	def __init__(self,max_ = 10,min_ = -10,maximal_nb_clusters = MAXNBCLUSTERS, maximal_nb_points_in_plot = MAXNBPOINTS, r_base=0.5, ncols=5, palette=[[int(255*c) for c in col] for col in sns.color_palette("hls", MAXNBCLUSTERS)]):
		PolarPlotter.__init__(self,max_,min_, only_medians=maximal_nb_points_in_plot == 0, r_base=r_base, ncols=ncols, palette=palette)
		self.plot = [pg.PlotItem(title="cluster"+str(_)) for _ in  xrange(maximal_nb_clusters)]
		self.maximal_distance = self.max-self.min+self.r_base
		for plot in self.plot:
			plot.setRange(xRange=[-self.maximal_distance,self.maximal_distance], yRange = [-self.maximal_distance,self.maximal_distance])
		self.data_points_plot = [pg.PlotDataItem(pen=pg.mkPen(color=(200, 200, 200,30))) for _ in xrange(maximal_nb_points_in_plot)]
		self.data_clusters_plot = [pg.PlotDataItem(pen=pg.mkPen(color=self.palette[i])) for i in xrange(maximal_nb_clusters)]
		self.plot_added = [0 for _ in  xrange(maximal_nb_clusters)]
		self.data = np.array([])
		for pl in self.plot:
			for point_plot in self.data_points_plot:
				 pl.addItem(point_plot)
		for i,cluster_plot in enumerate(self.data_clusters_plot):
			 self.plot[i].addItem(cluster_plot)	 
		self.oldest_point_index = 0
		self.maximal_nb_points_in_plot = maximal_nb_points_in_plot
		self.maximal_nb_clusters = maximal_nb_clusters
	
	def plot_dataset_points(self,dataset_points):
		for point in dataset_points:
			self.plot_point(point)
	
	def plot_point(self,point):
		point_plot_modified = -1
		if self.maximal_nb_points_in_plot > 0:
			if self.data.shape[0] < self.maximal_nb_points_in_plot or len(self.data.shape)==1:
				if self.data.shape[0] == 0 : 
					self.data = np.array(point)
					point_plot_modified = 0
				else:
					self.data = np.vstack((self.data,np.array(point)))
					point_plot_modified = self.data.shape[0]-1
			else:
				self.data[self.oldest_point_index] = np.array(point)
				point_plot_modified = self.oldest_point_index
				self.oldest_point_index += 1
				self.oldest_point_index = self.oldest_point_index % self.maximal_nb_points_in_plot
			self.n_dim = max(self.n_dim,len(point))
		#self.max = self.data.max().max()
		#self.min = self.data.min().min()
		x,y = self.point2shape(point)
		self.data_points_plot[point_plot_modified].setData(x,y)
	
	def plot_cluster(self,cluster_coordinates,cluster_id):
		x,y = self.point2shape(cluster_coordinates)
		self.data_clusters_plot[cluster_id].setData(x,y)
		self.win.addItem(self.plot[cluster_id])
		
	def plot_clusters(self,clusters_coordinates):
		for i,pl_added in enumerate(self.plot_added):
			if pl_added:
				self.win.removeItem(self.plot[i])
			self.plot_added[i] = 0
		for i,key in enumerate(clusters_coordinates.keys()):
			self.plot_cluster(clusters_coordinates[key],key)
			self.plot_added[key] = 1
			if len(clusters_coordinates.keys())/2 >0 and i%int(3*math.sqrt((len(clusters_coordinates.keys())/2))) == 0:
				self.win.nextRow()

if __name__ == '__main__':
    app = QtGui.QApplication([])

    path = 'wifi_transformed.h5'
    key = 'inv'

    df = pd.read_hdf(path, key)
    df.fillna(0,inplace=True)
    
    groups = range(10)
    df = df.loc[np.in1d(df.classes_fines,groups),:]
    
    data = df.iloc[:,18:-2]
    class_col = df.classes_fines

    plotter = PolarPlotter(data=data, 
                            class_col=class_col, 
                            only_medians=True,
                            )
    plotter.plot()

    QtGui.QApplication.instance().exec_()