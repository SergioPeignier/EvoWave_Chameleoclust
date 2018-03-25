# EvoWave Chameleoclust

![evowave visualization](./evowave_visualization)

## Chameleoclust

ChameleoClust is a bio-inspired algorithm implementing an evolvable genome structure, including several bio-like features such as a variable genome length, both functional and non-functional elements and mutation operators including chromosomal rearrangements.

The main purpose of the design of ChameleoClust+ is to take advantage of the large degree of freedom provided by its evolvable structure to detect various number of clusters in subspaces of various dimensions.
This algorithm has been developped in the context of the [EvoEvo project](https://evoevo.liris.cnrs.fr/)

## EvoWave

EvoWave, our first real-world application, deals with the Wi-Fi environment in which a micro-computer is immersed. This environment is defined by the strength of the sig- nal from every Wi-Fi antenna in the neighborhood. It depends especially on available routers and other computers or mobile devices, so it is linked to the context of use of the computer: work, teaching, house, etc. If the Wi-Fi signals from different contexts are dissimilar enough from each other, we expect that Chameleoclust should be able to discriminate these contexts using the data. This corresponds to a dynamic stream problem as new classes, i.e., context of use of the computers, are always susceptible to appear/disappear, and the present Wi-Fi antennas are never the same in different contexts (features appearing and disappearing). This example is also challenging regarding the high dimensionality and noise level of the data. The experimental setup is split in four steps:

+ Acquisition of data
+ Preprocessing
+ Clustering using Chameleoclust
+ Visualization of the resulting clusters and comparison with ground-truth


# INSTALLATION

Dependencies :
    - tshark (free and open-source. Official website : https://www.wireshark.org/download.html)
    - python 2.7 ([anaconda distribution by continuum analytics](https://www.continuum.io/downloads) contains every needed library, otherwise, you'll need to install numpy, pandas, seaborn and pyqtgraph)

Chameleoclust compilation :
    We provide binaries for both Linux and OS X. If you need to compile it yourself for any purpose, sources are provided, and can be compiled via the command "python2.7 setup.py build_ext --inplace" in the ChameleoClust/ folder.

# USAGE

+ `WifiCapture.py`
    Proceed to data acquisition. Full documentation is accessible through `python2.7 WifiCapture.py -h"

    Example : `python2.7 WifiCapture.py new my_name my_building my_floor my_office --colleague my_colleague --pathToLog ~/wifidata --capDuration 60 --capSpacing 120 --capNumber 10`
    This command will start a new acquisition identified by both the day and hour and info fields such as author, place, floor, room... Data will be stored in the `~/wifidata` folder. The acqusition will be split in ten parts of 1 minutes, spaced by 2-minutes pauses. Thus the total time of the acquisition will be 30 min

+ `WifiLog.py`
    Manage log files from `WifiCapture.py`. Full documentation is accessible through `python2.7 WifiLog.py -h`

    Conversion example from the project : `python WifiLog.py convert EVOWAVE_DATA/collect1 csv --override --macKeyLen 2`
    This command reads the raw log files and proceed to the dimensionality reduction step. The output is a csv file.

+ `csv2h5.py`
    Convert data from csv to pandas hdf5 format. Needed to use the scripts below.
    This is not a command line tool but a script. Input and output path are to be modified inside the script.

+ `preprocessing_logMax-log.py`
    Last steps of preprocessing.
    Input and output files are to be changed inside the script.

+ `evowave_clustering.py`
    Python script proceeding to the clustering step and building output needed for the visualization.
    Path to Chameleoclust lib, input file and parameters are to be changed inside the code. A simplified runable example of Chameleoclust+ is provided in the chameleo_example.py.

+ `video_evowave.py`
    Python script producing the visual from log files produced by the evowave_clustering.py script.

All scripts should be launched with python2.7
