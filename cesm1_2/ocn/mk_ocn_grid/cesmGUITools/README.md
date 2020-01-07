# cesmGUITools

cesmGUITools is a collection of python GUI applications meant to help edit geophysical fields used for the CESM model.

## Requirements
The python applications are written to make use of standard python packages as well as some widely used libraries. The following are python third-party libraries are the requirements for cesmGUITools:

1. Numpy 1.9.2
2. Matplotlib 1.4.3
3. Basemap 1.0.7
4. netcdf4-python 1.1.7
5. PyQt4

## Installation
There is no specific installation process required to use the scripts.

## Usage
All scripts are executable and meant to be called from the command-line with required and optional arguments. Calling help on any of the executable scripts on the command line by typing `./thisscript.py --help` will print the required and optional arguments that program accepts.

## Editors presently available
The following graphical editors are currently available:

1. TopoEditor

  This is used for editing a topography dataset.

2. KMTEditor

  This is used for editing the POP Ocean model's KMT levels.

## NOTE
If you get an error an error like this:
```
TypeError: 'PySide.QtGui.QWidget.setParent' called with wrong argument types:
  PySide.QtGui.QWidget.setParent(QWidget)
Supported signatures:
  PySide.QtGui.QWidget.setParent(PySide.QtGui.QWidget)
  PySide.QtGui.QWidget.setParent(PySide.QtGui.QWidget, PySide.QtCore.Qt.WindowFlags)
```
then it is because your ```matplotlibrc``` file is configured to use PySide backend for Qt4Agg. You have to change it to use the PyQt4 backend since that is the API that is being used in cesmGUITools editors. Basically you need to make sure that these two lines (uncommented) are in your ```matplotlibrc``` file:
```
backend      : Qt4Agg
backend.qt4 : PyQt4
```
See [this page](http://matplotlib.org/users/customizing.html) for information on where you can find ```matplotlibrc``` on your system.

acknowledgment:
Deepak Chandan, PhD Candidate
Earth Atmosphere and Planetary Physics Group
Department of Physics, University of Toronto
https://github.com/dchandan
