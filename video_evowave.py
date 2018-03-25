#-*- coding: utf-8 -*-

import itertools, time

import numpy as np
import pandas as pd
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore

class AutoPlotter(object):
    def __init__(self, frame_rate, polar_plotter, heat_map_plotter, info_collecter):
        self.frame_rate = frame_rate
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.update)

        self.polar_plotter = polar_plotter
        self.heat_map_plotter = heat_map_plotter
        self.info_collecter = info_collecter
        self.generation = self.info_collecter.generation

        print("début init_plot")
        self.view = pg.GraphicsView()
        self.l = pg.GraphicsLayout()
        self.view.setCentralItem(self.l)
        self.view.show()
        self.view.setWindowTitle('EvoWave : Supspace-clustering stream example')
        self.view.resize(1400,1400)
        self.build_plots_grid()
    
    def build_plots_grid(self):
        self.label = self.l.addLabel("Objects",size='15pt')
        self.l.nextRow()
        self.l.layout.setRowStretchFactor(1,1)
        self.l.layout.setRowStretchFactor(2,2)
        self.l.layout.setRowStretchFactor(3,1)

        self.l1 = self.l.addLayout(colspan=(self.polar_plotter.max_class/2+1), border=(50,0,0))
        self.l1.setContentsMargins(30, 0, 130, 0)

        for i in range(0,self.polar_plotter.max_class,2):
            self.polar_plotter.build_plot_nb(i, self.l1)

        self.l.nextRow()

        self.l2 = self.l.addLayout(colspan=(self.polar_plotter.max_class/2+1), border=(50,0,0))
        self.heat_map_plotter.build_plot(self.l2)
        self.l.nextRow()
        
        self.l3 = self.l.addLayout(colspan=(self.polar_plotter.max_class/2+1), border=(50,0,0))
        self.l3.setContentsMargins(160, 0, 0, 0)

        for i in range(1,self.polar_plotter.max_class,2):
            self.polar_plotter.build_plot_nb(i, self.l3)


    def update(self):
        self.update_data()
        self.update_plot()

    def update_data(self):
        info = self.info_collecter.info
        self.polar_plotter.update_data(model=info['model'],
                                        phenotype=info['phenotype'],
                                        )
        self.heat_map_plotter.update_data(confusion_matrix=info['confusion_matrix'])
        self.generation = info['generation']

    def update_plot(self):
        self.polar_plotter.update_plot()
        self.heat_map_plotter.update_plot()
        self.label.setText("Objects {} to {}".format(self.generation - 100, self.generation))

    def start(self):
        print(1000/self.frame_rate)
        self.timer.start(int(1000/self.frame_rate))
        print('okay')


class PlotBase(object):
    def __init__(self, win, max_, centroid, median, **kargs):
        self.plot = win.addPlot(**kargs)
        self.plot.setContentsMargins(45,15,45,15)
        self.plot.setRange(xRange=[-max_,max_], yRange = [-max_,max_])
        self.plot.enableAutoRange(False,False)
        self.plot.showAxis('left',False)
        self.plot.showAxis('bottom',False)

        self.curves = []
        self.centroid = centroid
        self.plot.addItem(self.centroid)
        self.median = median
        self.plot.addItem(self.median)

    def addItem(self, item):
        self.curves.append(item)
        self.plot.addItem(item)

    def clear(self):
        for curve in self.curves:
            self.plot.removeItem(curve)
        self.curves = []


class PolarPlotter(object):

    def __init__(self, data, class_col, r_base=0.5,  palette=list(itertools.permutations([50,150,200])) + list(itertools.permutations([0,100,200])), max_class=10):
        self.data = data
        self.class_col = class_col

        self.min = 0 
        self.max = self.data.max().max()

        self.max_class = max_class
        self.n_dim = self.data.shape[1] 
        self.classes = np.unique(self.class_col)

        self.r_base = r_base
        self.palette = palette
    
        self.plots = {}

    def polar2cartesian(self, r, dim):
        theta = (2 * np.pi * dim) / self.n_dim
        x = (r + self.r_base) * np.cos(theta)
        y = (r + self.r_base) * np.sin(theta)
        return x,y

    def point2shape(self, vec):
        return self.polar2cartesian(vec-self.min, np.arange(len(vec)))

    def draw_limit_circles(self, p):
        p.addItem(pg.PlotDataItem(x=[0,0],y=[0,0], pen=None, symbol='o', symbolBrush=None, symbolSize=[(self.min+self.r_base)*2, (self.max+self.r_base)*2], pxMode=False))


    def update_data(self, model, phenotype):
        self.data = model.iloc[:,:-2]
        self.class_col = model['cluster_found']
        self.classes = np.unique(self.class_col)
        #phenotype
        self.centroids = phenotype.fillna(0)
        self.medians = self.data.groupby(self.class_col).median()

    def update_plot(self):
        for i,plot in self.plots.items():
            plot.clear()
            # plot individual data points
            group = self.data.loc[self.class_col==i, :]#.iloc[::2]

            x, y = self.point2shape(group.median())
            plot.median.setPen(width = group.shape[0]*0.05)
            plot.median.setData(x,y)

        for idx, row in self.centroids.iterrows():
            try:
                self.plots[idx].centroid.setData(*self.point2shape(row))
            except KeyError, e:
                pass

    def build_plot_nb(self, class_, win):
        x, y  = self.point2shape(np.zeros(self.n_dim))
        self.plots[class_] = PlotBase(win, 
                                        max_=self.max + self.r_base, 
                                        title="Cluster {}".format(class_),
                                        centroid = pg.PlotDataItem(x, y, pen=pg.mkPen(color=self.palette[class_], width=3.5)),
                                        median = pg.PlotDataItem(x, y, pen=pg.mkPen(color=(200, 200, 200,255))),
                                        )
        self.draw_limit_circles(self.plots[class_].plot)



class HeatMapPlotter(object):
    def __init__(self, class_labels=None):
        pass

    def build_plot(self, win):
        self.plt = win.addPlot(title='Confusion matrix')
        self.img = pg.ImageItem()

        pos = np.array([0.0, 0.5, 1.0])
        color = np.array([[0,0,0,255], [128,20,20,255], [200,40,40,255]], dtype=np.ubyte)
        map = pg.ColorMap(pos, color)
        lut = map.getLookupTable(0.0, 1.0, 256)
        self.img.setLookupTable(lut)

        self.plt.addItem(self.img)
        self.plt.axes['left']['item'].setTicks(self.ticks())
        # self.plt.showGrid(True,True)
        self.build_grid()
        self.labels = {}
    # def update(self):
    #     self.update_data()
    #     self.update_plot()

    def build_grid(self):
        for h in range(10):
            self.horizontal_line(h)
        for h in [1,4,6,8,9]:
            self.horizontal_line(h,alpha=200)

    def horizontal_line(self, h, alpha=100):
        self.plt.addItem(pg.PlotDataItem(np.ones(11)*h, pen=pg.mkPen(color=(100, 100, 100,alpha))))

    def update_data(self, confusion_matrix):
        self.data = confusion_matrix

    def update_plot(self):
        self.img.setImage(self.data.values)

        self.clean_labels()
        max_y, max_x = self.data.shape
        for j in range(max_y):
            for i in range(max_x):
                if self.data.iloc[j,i]!=0:
                    try:
                        self.labels[(j,i)].setText(str(int(self.data.iloc[j,i])))
                    except KeyError:
                        label = pg.TextItem(str(int(self.data.iloc[j,i])), anchor=(-3,1.2))
                        self.plt.addItem(label)
                        label.setPos(j, i)
                        self.labels[(j,i)] = label

    def clean_labels(self):
        for label in self.labels.values():
            label.setText('')

    def ticks(self):
        ticks = [
                    [(10,"Outliers"), (9,"1"), (8,"2.A"), (6,"3.A"), (4,"4.A"), (1,"5")],
                    [(7,".B"), (5,".B"), (3,".B"), (2,".C")],
                ]
        return ticks


keygen = lambda seed, generation: "seed"+str(seed)+"/"+"generation"+str(generation)

class EvoWaveInfoCollecter(object):
    def __init__(self, model_file_name, phenotype_file_name, confusion_file_name, seed, generation_min):
        self.model_file_name = model_file_name
        self.phenotype_file_name =phenotype_file_name
        self.confusion_file_name = confusion_file_name
        self.seed = seed
        self.generation = generation_min
        self.confusion_base = pd.DataFrame(np.zeros(shape=(10,10)))
        self.confusion_base.columns = ['Outliers','1','2.A','2.B','3.A','3.B','4.A','4.B','4.C','5'][::-1]

    @property
    def info(self):
        self.key = keygen(self.seed,self.generation)
        info = {'model':self.model,'phenotype':self.phenotype,'confusion_matrix':self.confusion_matrix, 'generation':self.generation}
        self.generation += 1
        return info

    @property
    def model(self):
        return pd.read_hdf(self.model_file_name, self.key)

    @property
    def phenotype(self):
        return pd.read_hdf(self.phenotype_file_name, self.key)
    
    @property
    def confusion_matrix(self):
        self.update_confusion_base()
        confusion_matrix = pd.read_hdf(self.confusion_file_name, self.key)
        # print confusion_matrix
        self.confusion_base.update(confusion_matrix)
        return self.confusion_base

    def update_confusion_base(self):
        self.confusion_base.index = ["{}-{}".format(self.generation,i) for i in range(10)]
        self.confusion_base.iloc[:] = 0


if __name__ == '__main__':

 
    phenotype_file_name = "phenotype.h5"
    model_file_name = "model.h5"
    confusion_file_name = "confusion_matrix.h5"

    dataset = "datasets/wifi_filter_120_LogMax-LogX_final_classes.h5"
    precision = 10000

    df = pd.read_hdf(dataset)
    X = (df.iloc[:,-256-2:-2].fillna(0)*precision).astype(np.int)

    seed = 0
    generation_min = 101
 
    app = QtGui.QApplication([])

    info_collecter = EvoWaveInfoCollecter(model_file_name=model_file_name,
                                            phenotype_file_name=phenotype_file_name,
                                            confusion_file_name=confusion_file_name,
                                            seed=seed,
                                            generation_min=generation_min,
                                            )

    polar_plotter = PolarPlotter(data=X, 
                                    class_col=np.zeros(X.shape[0])-1, 
                                    r_base=20000,
                                    max_class=10,
                                    )

    heat_map_plotter = HeatMapPlotter()

    plotter = AutoPlotter(frame_rate=10,
                            polar_plotter=polar_plotter,
                            heat_map_plotter=heat_map_plotter,
                            info_collecter=info_collecter,
                            )

    time.sleep(16)
    plotter.start()

    app.exec_()   
