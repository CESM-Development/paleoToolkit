#!/usr/bin/env python

"""
RMaskEditor.py

This program allows one to edit a 2D KMT data pixel for pixel.

Author : Deepak Chandan
Date   : June 24th, 2015
"""

import sys, argparse, os
from netCDF4 import Dataset
import numpy as np
from PyQt4.QtCore import *
from PyQt4.QtGui import *


import matplotlib as mpl
from matplotlib import pylab as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

from PyCESM.utilities import nccopy

from matplotlib.widgets import Lasso
from matplotlib import path


mpl.rc('axes',edgecolor='w')


class LassoTool(object):
    def __init__(self, ax, nx, ny, callback):
        """
        ARGUMENTS
            nx  - number of columns in the view window
            ny  - number of rows in the view window
            callback - a function from the main program that the LassoTool class
                   calls when it has selected a list of points inside a lasso
        """
        self.axes   = ax
        self.canvas = ax.figure.canvas
        self.main_app_callback = callback
        
        yv, xv = np.mgrid[0:ny,0:nx]

        self.points = np.c_[xv.flatten(), yv.flatten()]
        del yv, xv
        # self.points contains coordinate locations of each pixel in the view
        # where the locations are indexed as (x_index, y_index), equivalently
        # as (column_index, row_index). 
        
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)


    def callback(self, verts):
        """
        This function is called by the matplotlib lasso widget when the lasso is released.
        ARGUMENTS
            verts - a list of tuples which collectively defines a path
        """
        p = path.Path(verts)  # We first create a closed path using the vertices
        
        # Selecting points that lie inside the lasso >>>>
        ind = p.contains_points(self.points)
        selected_points = self.points[ind,:]
        # <<<< Selecting points that lie inside the lasso
        
        self.canvas.draw_idle()
        # self.canvas.widgetlock.release(self.lasso)
        del self.lasso
        # Now that we are done with the lasso stuff, we call back an appropriate 
        # function in the main program (this function was passed as an argument
        # to the constructor for LassoTool)
        self.main_app_callback(selected_points)

    def onpress(self, event):
        """
        This function is called by matplotlib when a button is pressed on the view
        window.
        ARGUMENT
            event - matplotlib passes the event as the argument to this function
        """
        # if self.canvas.widgetlock.locked(): 
        #     print "here"
        #     return
        if event.inaxes is None: return
        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callback)
        # acquire a lock on the widget drawing
        # self.canvas.widgetlock(self.lasso)




class DataContainer(object):
    """
    DataContainer: A "data container" class for this application which does the job
    of storing the map and other related data as well as a "cursor" upon this data.
    """
    class Cursor(object):
        def __init__(self):
            self.marker = None   # The symbol which is used to denote the marker. Leave it as 'None' to show a square.
            self.x = 0           # X-position of the cursor
            self.y = 0           # Y-position of the cursor


    def __init__(self, nrows, ncols, fname, datavar):
        """
        ARGUMENTS
            nrows - number of rows for the view
            ncols - number of columns for the view
            fname - name of the data file.
            scale - a multiplicative scale factor for the data to be visualized and edited
        """
        self.fname   = fname
        self.datavar = datavar
        self.__read_nc_file()
        self.orig_data = np.copy(self.data)
        self.ny, self.nx = self.data.shape


        self.lon_modulo = 360

        # Datawindow variables
        self.view   = None        # The array (actually a numpy view) that stores the data to be displayed in the main window
        self.nrows  = nrows       # Number of rows to display in the main windows (the 'view')
        self.ncols  = ncols       # Number of cols to display in the main windows
        self.si     = None        # 0-based row index of the first element
        self.sj     = None        # 0-based col index of the first element

        # Tracking which elements are changed
        self.changes = np.zeros((self.ny*self.nx, 3), dtype=np.int32)
        self.changes_row_idx = 0

        # A cursor object on the view
        self.cursor = DataContainer.Cursor()


    def __read_nc_file(self):
        """ This subroutine reads the netCDF4 data file. """
        ncfile        = Dataset(self.fname, "r", format="NETCDF4")
        # We are flipping the arrays so that the latitudes go from 90:-90
        self.data     = np.flipud(ncfile.variables[self.datavar][:,:])
        self.kmt_lons = np.flipud(ncfile.variables["ULON"][:,:])
        self.kmt_lats = np.flipud(ncfile.variables["ULAT"][:,:])
        ncfile.close()

        # Now we are making the data array a masked array
        self.data = np.ma.array(self.data, mask=(self.data == 0))
        # Now we are making the mask of this masked array "hard", i.e. it cannot be changed
        self.data.harden_mask()
        self.data[:,:] = 50.0   # Putting some default value for the ocean regions



    def getViewStatistics(self): return (self.view.min(), self.view.max(), self.view.mean())


    def getCursor(self): return self.cursor

    def getValueUnderCursor(self):
        # We need to first map the position of the cursor in the view to the position of the
        # on the global map.
        ci, cj = self.viewIndex2GlobalIndex(self.cursor.y, self.cursor.x)
        return self.data[ci, cj]


    def updateCursorPosition(self, event):
        """
        Updates the current view in the data window. This function is
        triggered whenever the user presses the arrow keys.
        ARGUMENTS
            event - Qt key press event
        """
        key = event.key()
        if key == Qt.Key_Up:
            self.cursor.y = max(0, self.cursor.y - 1)
        elif key == Qt.Key_Down:
            self.cursor.y = min(self.nrows-1, self.cursor.y + 1)
        elif key == Qt.Key_Left:
            self.cursor.x = max(0, self.cursor.x - 1)
        elif key == Qt.Key_Right:
            self.cursor.x = min(self.ncols-1, self.cursor.x + 1)


    def updateView(self, si, sj):
        """
        Updates the data for the view.
        ARGUMENTS
            si, sj - the global 0-based i,j indices of the top left corner of the view
        RETURNS
            Statistics of the newly updated view (a tuple with the min, max, and the mean for the new view)
        """
        self.view = self.data[si:si+self.nrows, sj:sj+self.ncols].view()
        self.si   = si
        self.sj   = sj
        return self.getViewStatistics()



    def moveView(self, move):
        """
        Moves the view window over the global dataset in response to the L,R,U,D keys.
        ARGUMENTS
            move : A Qt key value
        RETURNS
            Statistics of the newly updated view (a tuple with the min, max, and the mean for the new view)
        """
        col_inc = int(self.ncols*0.25)  # Column increment
        row_inc = int(self.nrows*0.25)  # Row increment

        if (move == Qt.Key_L):  # MOVE THE VIEW RIGHT
            new_sj = min(self.nx-self.ncols, self.sj + col_inc)
            return self.updateView(self.si, new_sj)
        elif (move == Qt.Key_H): # MOVE THE VIEW LEFT
            new_sj = max(0, self.sj - col_inc)
            return self.updateView(self.si, new_sj)
        elif (move == Qt.Key_K): # MOVE THE VIEW UP
            new_si = max(0, self.si - row_inc)
            return self.updateView(new_si, self.sj)
        elif (move == Qt.Key_J): # MOVE THE VIEW DOWN
            new_si = min(self.ny-self.nrows, self.si + row_inc)
            return self.updateView(new_si, self.sj)


    def modifyValue(self, inp):
        """
        Modify the value for a particular pixel. The location of the pixel is that
        determined by the current position of the cursor.
        ARGUMENTS
            inp - a string containing a float
        """
        ci, cj = self.viewIndex2GlobalIndex(self.cursor.y, self.cursor.x)
        
        _tmp = int(float(inp))
        self.data[ci, cj] = _tmp
        self.changes[self.changes_row_idx, :] = ci, cj, _tmp
        self.changes_row_idx += 1



    def viewIndex2GlobalIndex(self, i, j):
        """ Converts an i,j index into the data window into an index for the
        same element into the global data. """
        return (self.si + i, self.sj + j)



class RMaskEditor(QMainWindow):

    def __init__(self, fname, datavar, dwx=60, dwy=60):
        """
        ARGUMENTS:
            fname    - Name of the netcdf4 file
            datavar  - Name of the data variable in the file for which to plot
            dwx, dwy - size of the DataContainer in number of array elements
            scale    - A float that will be multiplied with the data to scale the data
        """
        super(RMaskEditor, self).__init__(None)
        self.setWindowTitle('RMaskEditor - {0}'.format(fname))

        #  Creating a variable that contains all the data
        self.dc = DataContainer(dwy, dwx, fname, datavar)

        self.cursor = self.dc.getCursor()  # Defining a cursor on the data
        # This is the Rectangle boundary drawn on the world map that bounds the region
        # under the view. At the start of the program it's value is None
        self.prvrect = None

        # This stores the matplotlib text objects used in the rendering of the colorbar
        # This is used by the function draw_colorbar()
        self.colorbarlabels = []
        # self.colorbar_gradient is used by the function draw_colorbar()
        self.colorbar_gradient = np.vstack((np.linspace(0, 1, 256), np.linspace(0, 1, 256)))


        # The previously updated value
        self.buffer_value = None

        self.unsaved_changes_exist = False

        self.maps = mpl.cm.datad.keys()  # The names of colormaps available
        self.maps.sort() # Sorting them alphabetically for ease of use

        self.create_menu()
        self.create_main_window()

        self.set_stats_info(self.dc.updateView(0, 0))  # Set the initial view for the data container class

        self.draw_preview_worldmap()
        self.render_view()
        self.statusBar().showMessage('RMaskEditor 2015')


    def keyPressEvent(self, e):
        if e.key() == Qt.Key_Equal:
            # Pressing = for edit
            self.inputbox.setFocus()
        elif e.key() == Qt.Key_V:
            self.update_value(self.buffer_value)
        elif e.key() == Qt.Key_C:
            self.buffer_value = self.dc.getValueUnderCursor()
            self.statusBar().showMessage('Value copied: {0}'.format(self.buffer_value), 2000)
        elif e.key() == Qt.Key_A:
            # Get the 5-point average and update the value
            self.update_value_with_average(self.dc.getAverage())
        elif e.key() == Qt.Key_M:
            self.colormaps.setFocus()
        elif e.key() == Qt.Key_Escape:
            # Pressing escape to refocus back to the main frame
            self.main_frame.setFocus()
        elif e.key() in [Qt.Key_H, Qt.Key_J, Qt.Key_K, Qt.Key_L]:
            self.set_stats_info(self.dc.moveView(e.key()))
            self.render_view()
        else:
            self.dc.updateCursorPosition(e)
            self.draw_cursor()


    def create_main_window(self):
        """
        This function creates the main window of the program.
        """
        self.main_frame = QWidget()
        self.main_frame.setMinimumSize(QSize(700, 700))

        self.dpi = 100
        self.fig = plt.Figure((6, 6), dpi=self.dpi, facecolor='w', edgecolor='w')
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)


        # Stuff for the preview map >>>>>>>>>>>>>>>>>>>>>>>>>
        self.preview_frame = QWidget()
        self.preview_fig = plt.Figure((3, 1.6), dpi=self.dpi, facecolor='w', edgecolor='w')
        self.preview = FigureCanvas(self.preview_fig)
        self.preview.setParent(self.preview_frame)
        self.preview_axes = self.preview_fig.add_subplot(111)

        # self.preview_fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.preview_fig.subplots_adjust(top=1, bottom=0, left=0, right=1)
        # Stuff for the preview map <<<<<<<<<<<<<<<<<<<<<<<<<

        # Stuff for the colorbar >>>>>>>>>>>>>>>>>>>>>>>>>
        self.colorbar_frame = QWidget()
        self.colorbar_fig   = plt.Figure((3,0.4), dpi=self.dpi, facecolor='w', edgecolor='w')
        self.colorbar       = FigureCanvas(self.colorbar_fig)
        self.colorbar.setParent(self.colorbar_frame)
        self.colorbar_axes  = self.colorbar_fig.add_subplot(111)
        self.colorbar_fig.subplots_adjust(top=1.0, bottom=0.35, left=0.02, right=0.97) #Tightening the area around the subplot
        self.colorbar_fig.patch.set_facecolor('none') # Making the figure background transparent

        self.colorbar_axes.get_xaxis().set_visible(False)
        self.colorbar_axes.get_yaxis().set_visible(False)
        # Stuff for the colorbar <<<<<<<<<<<<<<<<<<<<<<<<<


        # Since we have only one plot, we can use add_axes
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)
        self.lman = LassoTool(self.axes, self.dc.ny, self.dc.nx, self.modify_selected_points)
        # Turning off the axes ticks to maximize space. Also the labels were meaningless
        # anyway because they were not representing the actual lat/lons.
        self.axes.get_xaxis().set_visible(False)
        self.axes.get_yaxis().set_visible(False)



        # Other GUI controls
        # Information section
        font = QFont("SansSerif", 14)

        # STATISTICS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        self.statdisplay = QLabel("Statistics:")
        self.statgrid    = QGridLayout()
        self.statgrid.setSpacing(5)
        w = QLabel("View"); w.setFont(font)
        self.statgrid.addWidget(w,  1, 1, Qt.AlignCenter)
        w = QLabel("Global"); w.setFont(font)
        self.statgrid.addWidget(w, 1, 2, Qt.AlignCenter)

        for i, name in enumerate(["Minimum", "Maximum", "Mean"]):
            w = QLabel(name)
            w.setFont(font)
            self.statgrid.addWidget(w, i+2, 0, Qt.AlignLeft)

        self.statsarray = []
        for i in range(6): self.statsarray.append(QLabel(''))

        self.statsarray[3].setText("{0:3d}".format(int(self.dc.data.min())))
        self.statsarray[4].setText("{0:3d}".format(int(self.dc.data.max())))
        self.statsarray[5].setText("{0:3d}".format(int(self.dc.data.mean())))

        for i in range(3):
            self.statgrid.addWidget(self.statsarray[i], i+2, 1, Qt.AlignCenter)
            self.statgrid.addWidget(self.statsarray[i+3], i+2, 2, Qt.AlignCenter)
            self.statsarray[i].setFont(font)
            self.statsarray[i+3].setFont(font)

        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        # PIXEL INFO >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        self.infodisplay = QLabel("Pixel Information:")
        self.infogrid    = QGridLayout()
        self.infogrid.setSpacing(5)


        rr = QLabel("")
        self.infogrid.addWidget(rr, 1, 1)

        for i, name in enumerate(["Lat", "Lon", "I, J", "Region"]):
            w = QLabel(name)
            w.setFont(font)
            # w.setLineWidth(1)
            # w.setFrameShape(QFrame.Box)
            self.infogrid.addWidget(w, i+2, 0, Qt.AlignLeft)

        self.latdisplay  = QLabel("")
        self.londisplay  = QLabel("")
        self.idxdisplay  = QLabel("")
        self.valdisplay  = QLabel("")
        for i,w in enumerate([self.latdisplay, self.londisplay, self.idxdisplay, self.valdisplay]):
            self.infogrid.addWidget(w, i+2, 1, Qt.AlignLeft)
            w.setFont(font)
        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        # Colorscheme selector
        cmap_label = QLabel('Colorscheme:')
        self.colormaps = QComboBox(self)
        self.colormaps.addItems(self.maps)
        self.colormaps.setCurrentIndex(self.maps.index('Set1'))
        self.colormaps.currentIndexChanged.connect(self.render_view)

        # New value editor
        valhbox = QHBoxLayout()
        w = QLabel("Enter new value: ")
        w.setFont(font)
        valhbox.addWidget(w)
        self.inputbox = QLineEdit()
        self.inputbox.returnPressed.connect(self.update_value)
        valhbox.addWidget(self.inputbox)

        for item in [self.statdisplay, self.infodisplay, self.latdisplay, \
                     self.idxdisplay, self.londisplay, \
                     self.valdisplay, cmap_label]:
            item.setFont(font)


        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)


        vbox2 = QVBoxLayout()
        vbox2.addWidget(self.preview)           # Adding preview window

        # Horizontal line
        fr = QFrame(); fr.setLineWidth(0.5); fr.setFrameShape(QFrame.HLine)
        vbox2.addWidget(fr)
        vbox2.setAlignment(fr, Qt.AlignTop)


        # databox is a horizontal box for the displayed "information"
        databox  = QHBoxLayout()
        statvbox = QVBoxLayout()    # This is the sub-box for the statistics
        statvbox.addWidget(self.statdisplay)
        statvbox.setAlignment(self.statdisplay, Qt.AlignTop)
        statvbox.addLayout(self.statgrid)
        databox.addLayout(statvbox)

        fr = QFrame(); fr.setLineWidth(0.5); fr.setFrameShape(QFrame.VLine)
        databox.addWidget(fr)

        pixelvbox = QVBoxLayout()   # This is another sub-box of databox for the pixel info
        pixelvbox.addWidget(self.infodisplay)
        pixelvbox.setAlignment(self.infodisplay, Qt.AlignTop)
        pixelvbox.addLayout(self.infogrid)
        databox.addLayout(pixelvbox)

        vbox2.addLayout(databox)

        # Horizontal line
        fr = QFrame(); fr.setLineWidth(0.5); fr.setFrameShape(QFrame.HLine)
        vbox2.addWidget(fr)
        vbox2.setAlignment(fr, Qt.AlignTop)



        vbox2.addLayout(valhbox, Qt.AlignTop)

        vbox2.addStretch(1)
        vbox2.addWidget(cmap_label)           # Adding the colormap label
        vbox2.addWidget(self.colormaps)       # Adidng the colormap selector
        vbox2.addWidget(self.colorbar)
        vbox2.setAlignment(self.colorbar, Qt.AlignCenter)
        vbox2.addStretch(1)

        helpbox = QVBoxLayout()

        help_title = QLabel("Help")
        font = QFont("SansSerif", 14)
        font.setBold(True)
        help_title.setFont(font)

        helpbox.addWidget(help_title); helpbox.setAlignment(help_title, Qt.AlignTop)

        helpgrid = QGridLayout()
        helpgrid.setSpacing(5)

        helpgrid.addWidget(QLabel("Key"),                 0, 0, 1, 1,Qt.AlignLeft)
        helpgrid.addWidget(QLabel("Action"),              0, 1, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("up, down,left,right"), 1, 0, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("move cursor"),         1, 1, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("h,j,k,l"),             2, 0, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("move view"),           2, 1, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("a"),                   3, 0, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("replace cell value with 9 point average"),         3, 1, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("="),                   4, 0, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("enter new value for cell"),         4, 1, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("c"),                   5, 0, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("copy value under cursor into buffer"),         5, 1, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("v"),                   6, 0, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("paste value in buffer"),6, 1, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("m"),                   7, 0, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("move focus to color selector"),         7, 1, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("Escape"), 8, 0, 1, 1, Qt.AlignLeft)
        helpgrid.addWidget(QLabel("move focus to main view"),         8, 1, 1, 1, Qt.AlignLeft)



        helpbox.addLayout(helpgrid); helpbox.setAlignment(helpgrid, Qt.AlignTop)


        vbox2.addLayout(helpbox)
        vbox2.setAlignment(helpbox, Qt.AlignBottom)


        hbox = QHBoxLayout()
        hbox.addLayout(vbox)
        hbox.addLayout(vbox2)

        self.main_frame.setLayout(hbox)
        self.setCentralWidget(self.main_frame)
        self.main_frame.setFocus()


    def draw_preview_worldmap(self):
        """
        This function draws the world map in the preview window on the top right hand corner
        of the application.
        """
        # cmap = mpl.cm.get_cmap(self.maps[self.colormaps.currentIndex()])
        self.preview_axes.imshow(self.dc.data, cmap=mpl.cm.Dark2, vmin=0.0, vmax=50.0, interpolation="none", extent=[0, 360, 0, 180])
        self.preview.draw()


    def draw_colorbar(self):
        """
        This function draws the colorbar and the labels for the colorbar. 
        """
        # First clear the colorbar axes
        self.colorbar_axes.cla()
        # Get the current selected colormap
        cmap = mpl.cm.get_cmap(self.maps[self.colormaps.currentIndex()])
        # Plot the colormap
        self.colorbar_axes.imshow(self.colorbar_gradient, aspect='auto', cmap=cmap)
        pos = list(self.colorbar_axes.get_position().bounds)
        dx = (pos[2])/6.  # The spacing between the labels

        if self.colorbarlabels != []:
            # Delete and clear the text labels
            while self.colorbarlabels: self.colorbarlabels.pop().remove()
        
        for i,l in enumerate([1,10,20,30,40,50,60]):
            self.colorbarlabels.append(self.colorbar_fig.text(pos[0]+(i*dx), 0.01, str(l), va='bottom', ha='center', fontsize=10))
        self.colorbar.draw()



    def draw_cursor(self, noremove=False):
        if self.cursor.marker and (not noremove): self.cursor.marker.remove()
        # The increment by 0.5 below is done so that the center of the marker is shifted
        # so that the top left corner of the square marker conicides with the top right corner
        # of each pixel. The value of 0.5 comes simply because each pixel are offset by 1 in each dimension.
        _cx, _cy = self.cursor.x+0.5, self.cursor.y+0.5
        self.cursor.marker = self.axes.scatter(_cx, _cy, s=55,
                             marker='s', edgecolor="k", facecolor='none', linewidth=2)
        self.set_information(self.cursor.y, self.cursor.x)
        self.canvas.draw()


    def render_view(self, clear=True):
        """
        ARGUMENTS
            clear  -  If true, then the axes clear method is called, else it is not called
        """
        self.draw_colorbar()
        if clear: self.axes.clear()
        # Either select the colormap through the combo box or specify a custom colormap
        # cmap = mpl.cm.get_cmap(self.maps[self.colormaps.currentIndex()])
        self.axes.pcolor(self.dc.view, cmap=mpl.cm.Dark2, edgecolors='k', linewidths=0.5, vmin=0.0, vmax=50.0)

        tmp1 = self.dc.nrows
        tmp2 = self.dc.ncols

        # Setting the axes limits. This helps in setting the right orientation of the plot
        # and in clontrolling how much extra space we want around the scatter plot.
        # I am putting 4% space around the scatter plot
        self.axes.set_ylim([int(tmp1*1.02), 0 - int(tmp1*0.02)])
        self.axes.set_xlim([0 - int(tmp2*0.02), int(tmp2*1.02)])
        self.canvas.draw()
        self.fig.tight_layout()
        self.draw_cursor(noremove=clear)




    def update_value(self, inp=None):
        if inp == None:
            inp = self.inputbox.text()   # Get the value in the text box
            if self.InputValueIsAcceptable(inp):
                self.dc.modifyValue(inp)     # Modify the data array
                self.unsaved_changes_exist = True
                self.buffer_value = inp
                self.statusBar().showMessage('Value changed: {0}'.format(inp), 2000)
                self.set_stats_info(self.dc.getViewStatistics())
                self.inputbox.clear()        # Now clear the input box
                self.render_view()           # Render the new view (which now contains the updated value)
                self.main_frame.setFocus()   # Bring focus back to the view



    def InputValueIsAcceptable(self, inp):
        """ This functions checks to see if the user has entered a valid
            vlaue kmt value. Returns true if the input value is acceptable. """
        tt = float(inp)
        if (tt != int(tt)):
            msg = "KMT level must be an integer value\nYou entered {0}".format(inp)
            QMessageBox.critical(self, "Invalid Value", msg.strip())
            return False

        if ((tt < 0) or (tt > 60)):
            msg = "KMT level values must be between 0 and 60.\nYou entered {0}".format(inp)
            QMessageBox.critical(self, "Invalid Value", msg.strip())
            return False
        return True



    def set_information(self, i, j):
        """ Sets the displayed information about the pixel in the right sidebar.
        ARGUMENTS
            i, j : the local (i.e. DataContainer) 0-based indices for the element
        """
        i_global, j_global = self.dc.viewIndex2GlobalIndex(i, j) # Convert local indices to global indices
        self.latdisplay.setText("{0:7.3f}".format(self.dc.kmt_lats[i_global, j_global]))
        self.londisplay.setText("{0:7.3f}".format(self.dc.kmt_lons[i_global, j_global]))
        self.idxdisplay.setText("{0:3d},{1:3d}".format(int(i_global), int(j_global)))
        try:
            # This will work if the cursor is not over a masked region. Otherwise the returned
            # value will be a masked value and the conversion to int will throw an exception
            self.valdisplay.setText("{0:3d}".format(int(self.dc.data[i_global, j_global])))
        except np.ma.core.MaskError:
            self.valdisplay.setText("{0:3d}".format(int(self.dc.data[i_global, j_global].item())))


    def set_stats_info(self, s):
        """
        Updates the statistics display panel with the stats for the view.
        ARGUMENTS
            s - a tuple with the min, max and mean of the view
        """
        self.statsarray[0].setText("{0:3d}".format(int(s[0])))
        self.statsarray[1].setText("{0:3d}".format(int(s[1])))
        self.statsarray[2].setText("{0:3d}".format(int(s[2])))



    def on_about(self):
        msg = """ Edit 2D geophysical field.  """
        QMessageBox.about(self, "About", msg.strip())


    def save_data(self):
        """
        Saves the data to the netCDF4 file from which input data was read in. First, the program asks
        for a variable name to use, but it remembers this variable name for subsequent saves.
        """

        # We construct a standard output name of the form: inputname_fixed.nc
        ofile = self.dc.fname.split(".")[0]+"_fixed.nc"

        clobber = False

        # Now we check if the file exists. Then we must ask the user if she wants to overwrite the file
        if os.path.exists(ofile):
            reply = QMessageBox.question(self, 'Message', "An output file exists. Do you want to overwrite?",
                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

            if reply == QMessageBox.Yes:
                clobber=True
            else:
                clobber=False

        if not clobber:
            # The use has responded to not overwrite, we must now ask for a new filename
            ofile, ok = QInputDialog.getText(self, "Saving KMT edit", "Enter name of a different output file:",)
            ofile = str(ofile)
            if (not ok):
                self.statusBar().showMessage('Save cancelled', 2000)
                return

        # If all is well so far, we have an acceptable output filename and we can proceed with writing.

        # To simplify creation of the output filename, i first duplicate the original file, then overwrite
        # the kmt array.
        nccopy(self.dc.fname, ofile, quiet=True, clobber=True, zlib=False, shuffle=False, classic=1)

        ncfile = Dataset(ofile, "a", format="NETCDF4")
        ncfile.variables["kmt"][:,:] = np.flipud(self.dc.data)
        self.unsaved_changes_exist = False
        self.statusBar().showMessage('Saved to file: %s' % ofile, 2000)




    def create_action(self, text, slot=None, shortcut=None,
                      icon=None, tip=None, checkable=False,
                      signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action


    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)


    def create_menu(self):
        self.file_menu = self.menuBar().addMenu("&File")

        load_file_action = self.create_action("&Save Data",
            shortcut="Ctrl+S", slot=self.save_data,
            tip="Save the data array")

        self.add_actions(self.file_menu, (load_file_action,))

        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About",
            shortcut='F1', slot=self.on_about,
            tip='About RMaskEditor')
        self.add_actions(self.help_menu, (about_action,))


    def closeEvent(self, event):
        if self.unsaved_changes_exist:
            reply = QMessageBox.question(self, 'Message', "There are unsaved changes. Are you sure you want to quit?",
                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

            if reply == QMessageBox.Yes:
                event.accept()
            else:
                event.ignore()


    def modify_selected_points(self, points):
        """
        ARGUMENTS
            points - 
        """
        # First we convert from view coordinates to global coordinates
        points_i, points_j = self.dc.viewIndex2GlobalIndex(points[:,1], points[:,0])

        # We prompt user to enter a value to define the oceanic region
        val, ok = QInputDialog.getText(self, "", "Enter oceanic region value:",)
        if (not ok):
            self.statusBar().showMessage('No ocean cells changed', 2000)
            return
        else:
            self.dc.data[points_i,points_j] = int(str(val))
            self.render_view(clear=False)

            self.draw_preview_worldmap()  # We update the preview map
            self.statusBar().showMessage('{0} ocean cells changed'.format(len(points_i), 2000))




def main():
    app = QApplication([])   # Create an application
    app.setWindowIcon(QIcon('Resources/rmaskicon.png'))

    parser = argparse.ArgumentParser(description='RMaskEditor', add_help=False)
    parser.add_argument('fname', nargs=1, type=str, help='name of the netcdf4 data file')
    parser.add_argument('-s',    nargs=1, type=int, help='size of the view in number of pixels', default=[60])
    args = parser.parse_args()

    mw = RMaskEditor(args.fname[0], "kmt", dwx=args.s[0], dwy=args.s[0])
    mw.show()     # Render the window
    mw.raise_()   # Bring the PyQt4 window to the front
    app.exec_()   # Run the application loop


if __name__ == "__main__":
    main()
