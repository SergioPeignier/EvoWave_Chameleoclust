# INSTALLATION

Dependencies :
    - tshark (free and open-source. Official website : https://www.wireshark.org/download.html)
    - python 2.7 (anaconda distribution by continuum analytics (https://www.continuum.io/downloads) contains every needed library, otherwise, you'll need to install numpy, pandas, seaborn and pyqtgraph)

Chameleoclust compilation :
    We provide binaries for both Linux and OS X. If you need to compile it yourself for any purpose, sources are provided, and can be compiled via the command "python2.7 setup.py build_ext --inplace" in the ChameleoClust/ folder.

# USAGE

    - WifiCapture.py
        Proceed to data acquisition. Full documentation is accessible through "python2.7 WifiCapture.py -h"

        Example : python2.7 WifiCapture.py new my_name my_building my_floor my_office --colleague my_colleague --pathToLog ~/wifidata --capDuration 60 --capSpacing 120 --capNumber 10
        This command will start a new acquisition identified by both the day and hour and info fields such as author, place, floor, room... Data will be stored in the ~/wifidata folder. The acqusition will be split in ten parts of 1 minutes, spaced by 2-minutes pauses. Thus the total time of the acquisition will be 30 min

    - WifiLog.py
        Manage log files from WifiCapture.py. Full documentation is accessible through "python2.7 WifiLog.py -h"

        Conversion example from the project : python WifiLog.py convert EVOWAVE_DATA/collect1 csv --override --macKeyLen 2
        This command reads the raw log files and proceed to the dimensionality reduction step. The output is a csv file.

    - csv2h5.py
        Convert data from csv to pandas hdf5 format. Needed to use the scripts below.
        This is not a command line tool but a script. Input and output path are to be modified inside the script.

    - preprocessing_logMax-log.py
        Last steps of preprocessing.
        Input and output files are to be changed inside the script.

    - evowave_clustering.py
        Python script proceeding to the clustering step and building output needed for the visualization.
        Path to Chameleoclust lib, input file and parameters are to be changed inside the code. A simplified runable example of Chameleoclust+ is provided in the chameleo_example.py.

    - video_evowave.py
        Python script producing the visual from log files produced by the evowave_clustering.py script.

All scripts should be launched with python2.7
