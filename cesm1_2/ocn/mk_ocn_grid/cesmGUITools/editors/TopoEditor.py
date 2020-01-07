#!/usr/bin/env python

"""
TopoEditor.py

This program allows one to edit a 2D latitude-longitude gridded
data pixel for pixel.

Author : Deepak Chandan
Date   : February 17th, 2015
"""

import sys, argparse
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

from cesmGUITools.utilities.topoutils import topography_cmap, make_balanced

mpl.rc('axes',edgecolor='w')

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


	def __init__(self, nrows, ncols, fname, datavar, scale):
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

		self.scale = scale
		self.data*=scale
		
		# Determining whether the longitude ranges from -180 to 180 or 0 to 360
		# this will determine how we plot the preview plot
		if (self.lons.min() < 0) and (self.lons.max() > 0):
			self.lon_modulo = 180
		else:
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
		""" This subroutine reads the netCDF4 data file. It looks for common names
		of the latitude and longitude variables in the file. If it cannot find any
		one of these coordinates, then it raises and error. """
		self.lons = None
		self.lats = None
		ncfile = Dataset(self.fname, "r", format="NETCDF4")
		for var in ["longitudes", "longitude", "lons"]:
			try:
				self.lons = ncfile.variables[var][:]
				self.lon_var = ncfile.variables[var].dimensions[0]
			except:
				pass

		for var in ["latitudes", "latitude", "lats"]:
			try:
				self.lats = ncfile.variables[var][:]
				self.lat_var = ncfile.variables[var].dimensions[0]
			except:
				pass

		if (self.lats == None) or (self.lons == None):
			ncfile.close()
			QMessageBox.critical(QWidget(), 'Error', "Latitude and or longitude variables not found.", QMessageBox.Ok)
			sys.exit()
			
		self.data = ncfile.variables[self.datavar][:,:]
		ncfile.close()



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
		# on_boundary = None   # This will store which boundary if any we've reached
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
		self.si = si
		self.sj = sj
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

	
	def modifyValue(self, input):
		"""
		Modify the value for a particular pixel. The location of the pixel is that
		determined by the current position of the cursor.
		ARGUMENTS
			input - a string containing a float (note: no data validity check is 
					performed on this string)
		"""
		ci, cj = self.viewIndex2GlobalIndex(self.cursor.y, self.cursor.x)

		_tmp = float(input)
		self.data[ci, cj] = _tmp
		self.changes[self.changes_row_idx, :] = ci, cj, _tmp
		self.changes_row_idx += 1



	def getAverage(self, center=False):
		"""
		Returns the average value at the cursor computed from the values of the surrounding cells. This
		is a 8 point average.
		If center==True, then includes the cell at which the cursor is in the calculation of the average,
		in which case it will be a 9 point average. 
		"""
		ci, cj = self.viewIndex2GlobalIndex(self.cursor.y, self.cursor.x)
		_sum   = self.data[ci-1:ci+2,cj-1:cj+2].sum()
		if center: 
			_sum += self.data[ci,cj]
			return _sum/9.
		else:
			_sum -= self.data[ci,cj]
			return _sum/8.


	def viewIndex2GlobalIndex(self, i, j):
		""" Converts an i,j index into the data window into an index for the
		same element into the global data. """
		return (self.si + i, self.sj + j)
	
	

class TopoEditor(QMainWindow):    

	def __init__(self, fname, datavar, dwx=60, dwy=60, scale=1.0):
		"""
		ARGUMENTS:
			fname    - Name of the netcdf4 file
			datavar  - Name of the data variable in the file for which to plot
			dwx, dwy - size of the DataContainer in number of array elements
			scale    - A float that will be multiplied with the data to scale the data
		"""
		super(TopoEditor, self).__init__(None)
		self.setWindowTitle('TopoEditor - {0}'.format(fname))
		
		#  Creating a variable that contains all the data
		self.dc = DataContainer(dwy, dwx, fname, datavar, scale)
		
		self.cursor = self.dc.getCursor()  # Defining a cursor on the data
		# This is the Rectangle boundary drawn on the world map that bounds the region
		# under the view. At the start of the program it's value is None
		self.prvrect = None                

		# The previously updated value
		self.buffer_value = None

		self.unsaved_changes_exist = False

		# The netcdf variable name for saving the modified data. The user will be asked
		# to enter the value when saving. 
		self.save_var   = None

		self.maps = mpl.cm.datad.keys()  # The names of colormaps available
		self.maps.sort() # Sorting them alphabetically for ease of use

		self.create_menu()
		self.create_main_window()

		self.set_stats_info(self.dc.updateView(0, 0))  # Set the initial view for the data container class

		self.draw_preview_worldmap()
		self.render_view()
		self.statusBar().showMessage('TopoEditor 2015')
	
	
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
			self.update_value_2(self.dc.getAverage(center=True))
		elif e.key() == Qt.Key_F:
			# Get the 4-point average and update the value
			self.update_value_2(self.dc.getAverage())
		# elif e.key() == Qt.Key_C:
		#     self.colormaps.setFocus()
		elif e.key() == Qt.Key_Escape:
			# Pressing escape to refocus back to the main frame
			self.main_frame.setFocus()
		elif e.key() in [Qt.Key_H, Qt.Key_J, Qt.Key_K, Qt.Key_L]:
			self.set_stats_info(self.dc.moveView(e.key()))
			self.render_view()
			self.render_edited_cells()
			self.draw_preview_rectangle()
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
		
		
		self.preview_frame = QWidget()
		self.preview_fig = plt.Figure((3, 1.6), dpi=self.dpi, facecolor='w', edgecolor='w')
		self.preview = FigureCanvas(self.preview_fig)
		self.preview.setParent(self.preview_frame)
		self.preview_axes = self.preview_fig.add_subplot(111)
		
		self.preview_fig.canvas.mpl_connect('button_press_event', self.onclick)
		self.preview_fig.subplots_adjust(top=1, bottom=0, left=0, right=1)
		
		# Since we have only one plot, we can use add_axes 
		# instead of add_subplot, but then the subplot
		# configuration tool in the navigation toolbar wouldn't
		# work.
		#
		self.axes = self.fig.add_subplot(111)
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
		w = QLabel("Local"); w.setFont(font)
		self.statgrid.addWidget(w,  1, 1, Qt.AlignCenter)
		w = QLabel("Global"); w.setFont(font)
		self.statgrid.addWidget(w, 1, 2, Qt.AlignCenter)

		for i, name in enumerate(["Minimum", "Maximum", "Mean"]):
			w = QLabel(name)
			w.setFont(font)
			self.statgrid.addWidget(w, i+2, 0, Qt.AlignLeft)

		self.statsarray = []
		for i in range(6): self.statsarray.append(QLabel(''))

		self.statsarray[3].setText("{0:5.2f}".format(self.dc.data.min()))
		self.statsarray[4].setText("{0:5.2f}".format(self.dc.data.max()))
		self.statsarray[5].setText("{0:5.2f}".format(self.dc.data.mean()))

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

		for i, name in enumerate(["Latitude", "Longitude", "Value"]):
			w = QLabel(name)
			w.setFont(font)
			self.infogrid.addWidget(w, i+1, 0, Qt.AlignLeft)

		self.latdisplay  = QLabel("")
		self.londisplay  = QLabel("")
		self.valdisplay  = QLabel("")
		for i,w in enumerate([self.latdisplay, self.londisplay, self.valdisplay]):
			self.infogrid.addWidget(w, i+1, 1, Qt.AlignLeft)
		# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		
		# Colorscheme selector
		cmap_label = QLabel('Colorscheme:')
		self.colormaps = QComboBox(self)
		self.colormaps.addItems(self.maps)
		self.colormaps.setCurrentIndex(self.maps.index('Spectral'))
		self.colormaps.currentIndexChanged.connect(self.render_view)

		# New value editor
		hbox = QHBoxLayout()
		w = QLabel("Enter new value: ")
		w.setFont(font)
		hbox.addWidget(w)
		self.inputbox = QLineEdit()
		self.inputbox.returnPressed.connect(self.update_value)
		hbox.addWidget(self.inputbox)

		for item in [self.statdisplay, self.infodisplay, self.latdisplay, self.londisplay, self.valdisplay, cmap_label]:
			item.setFont(font)

		
		vbox = QVBoxLayout()
		vbox.addWidget(self.canvas)
		
		
		vbox2 = QVBoxLayout()
		vbox2.addWidget(self.preview)
		
		vbox2.addWidget(self.statdisplay)
		vbox2.setAlignment(self.statdisplay, Qt.AlignTop)
		vbox2.addLayout(self.statgrid)

		vbox2.addWidget(self.infodisplay)
		vbox2.setAlignment(self.infodisplay, Qt.AlignTop)
		vbox2.addLayout(self.infogrid)

		# vbox2.addWidget(self.inputbox, Qt.AlignTop)
		vbox2.addLayout(hbox, Qt.AlignTop)

		vbox2.addStretch(1)
		vbox2.addWidget(cmap_label)
		vbox2.addWidget(self.colormaps)
		vbox2.addStretch(1)
			
		
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
		if self.dc.lon_modulo == 180:
			m = Basemap(projection='cyl', lon_0=0, llcrnrlat=-90,urcrnrlat=90,\
				llcrnrlon=-180,urcrnrlon=180,resolution='c', ax=self.preview_axes)
			self.preview_axes.set_xlim([-180,180])
		else:
			m = Basemap(projection='cyl', lon_0=180, llcrnrlat=-90,urcrnrlat=90,\
				llcrnrlon=0,urcrnrlon=360,resolution='c', ax=self.preview_axes)
			self.preview_axes.set_xlim([0,360])
		
		m.drawcoastlines(linewidth=0.5)
		m.fillcontinents()
		self.preview_axes.set_ylim([-90,90])
		self.draw_preview_rectangle()
	
	
	def draw_preview_rectangle(self):
		"""
		This function draws the Rectangle, which indicates the current region being shown
		in the view, in the preview window.
		"""
		# If a rectangle object exists, then remove it before drawing a new one
		if self.prvrect: self.prvrect.remove()

		# rect_llc_x and rect_llc_y are the x and y values of the lower left corner of the preview rectangle.
		rect_llc_x = self.dc.lons[self.dc.sj]
		rect_llc_y = self.dc.lats[min(self.dc.si+self.dc.nrows-1, self.dc.ny-1)]

		# dlon and dlat are the width and height of the preview rectangle
		dlon = abs(rect_llc_x - self.dc.lons[min(self.dc.sj+self.dc.ncols-1, self.dc.nx-1)])
		dlat = self.dc.lats[self.dc.si] - rect_llc_y
		
		self.prvrect = mpatches.Rectangle((rect_llc_x, rect_llc_y), dlon, dlat, linewidth=1, facecolor='g', alpha=0.3)
		self.preview_axes.add_patch(self.prvrect)
		self.preview.draw()
		
	
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
	

	def render_edited_cells(self):
		"""
		This function draws a box around cells that have been edited, thereby highlighting
		them from the other cells. 
		"""
		changes_row_idx = self.dc.changes_row_idx

		# We only need to go ahead if a change has been made, i.e. changes_row_idx is larger than 0
		if changes_row_idx > 0:
			# We are creating simple views to make the code for this function simpler
			change_rows = self.dc.changes[:self.dc.changes_row_idx, 0].view()  # A view on rows indices that have been changed
			change_cols = self.dc.changes[:self.dc.changes_row_idx, 1].view()  # A view on column indices that have been changed

			# Retrieving some class variables and storing for speed
			si    = self.dc.si
			sj    = self.dc.sj
			nrows = self.dc.nrows
			ncols = self.dc.ncols

			# We select the indices that lie within the view box
			indices_of_interest = np.logical_and(np.logical_and(change_rows >= si, 
																change_rows <= (si+nrows)), 
												 np.logical_and(change_cols >= sj, 
																change_cols <= (sj+ncols)))

			_i = change_rows[indices_of_interest] + 0.5 - si
			_j = change_cols[indices_of_interest] + 0.5 - sj
			
			self.axes.scatter(_j, _i, s=38, marker='s', edgecolor="k", facecolor='none', linewidth=1)
			self.canvas.draw()

	
	
	def render_view(self):
		self.axes.clear()
		# Either select the colormap through the combo box or specify a custom colormap
		# cmap = mpl.cm.get_cmap(self.maps[self.colormaps.currentIndex()])
		cmap   = topography_cmap(80, end=0.85)
		ll, ul = make_balanced(ll=-7.*self.dc.scale)
		self.axes.pcolor(self.dc.view, cmap=cmap, edgecolors='w', linewidths=0.5, vmin=ll, vmax=ul)
		
		# Setting the axes limits. This helps in setting the right orientation of the plot
		# and in clontrolling how much extra space we want around the scatter plot.
		tmp1 = self.dc.nrows
		tmp2 = self.dc.ncols
		# I am putting 4% space around the scatter plot
		self.axes.set_ylim([int(tmp1*1.02), 0 - int(tmp1*0.02)])
		self.axes.set_xlim([0 - int(tmp2*0.02), int(tmp2*1.02)])
		self.canvas.draw()
		self.fig.tight_layout()
		self.draw_cursor(noremove=True)
	
	


	def update_value(self, inp=None):
		if inp == None:
			inp = self.inputbox.text()   # Get the value in the text box
		self.dc.modifyValue(inp)     # Modify the data array
		self.unsaved_changes_exist = True 
		self.buffer_value = inp
		self.statusBar().showMessage('Value changed: {0}'.format(inp), 2000)
		self.set_stats_info(self.dc.getViewStatistics()) 
		self.inputbox.clear()        # Now clear the input box
		self.render_view()           # Render the new view (which now contains the updated value)
		self.render_edited_cells()
		self.main_frame.setFocus()   # Bring focus back to the view
	
	
	def update_value_2(self, val):
		"""
		Updates the value at the current cursor position using the input val.
		"""
		self.dc.modifyValue(val)     # Modify the data array
		self.unsaved_changes_exist = True 
		self.statusBar().showMessage('Value changed: {0}'.format(val), 2000)
		self.set_stats_info(self.dc.getViewStatistics()) 
		self.render_view()           # Render the new view (which now contains the updated value)   
		self.render_edited_cells()         

		
	
	def set_information(self, i, j):
		""" Sets the displayed information about the pixel in the right sidebar. 
		ARGUMENTS
			i, j : the local (i.e. DataContainer) 0-based indices for the element
		"""
		i_global, j_global = self.dc.viewIndex2GlobalIndex(i, j) # Convert local indices to global indices
		self.latdisplay.setText("{0}".format(self.dc.lats[i_global]))
		self.londisplay.setText("{0}".format(self.dc.lons[j_global]))
		self.valdisplay.setText("{0}".format(self.dc.data[i_global, j_global]))
	

	def set_stats_info(self, s):
		"""
		Updates the statistics display panel with the stats for the view.
		ARGUMENTS
			s - a tuple with the min, max and mean of the view
		"""
		self.statsarray[0].setText("{0:5.2f}".format(s[0]))
		self.statsarray[1].setText("{0:5.2f}".format(s[1]))
		self.statsarray[2].setText("{0:5.2f}".format(s[2]))


	def onclick(self, event):
		# 1. Get the global row, column indices of the point where mouse was clicked
		if (event.xdata == None) or (event.ydata == None): return
		px = np.where(abs(self.dc.lons - event.xdata) < 0.5)[0][0]
		py = np.where(abs(self.dc.lats - event.ydata) < 0.5)[0][0]
		# 2. Update the view data array 
		self.set_stats_info(self.dc.updateView(py, px))
		# 3. Render the view
		self.render_view()
		# 4. Set the cursor to be at the top-left corner
		self.cursor.x = 0
		self.cursor.y = 0
		# 5. Draw the cursor
		self.draw_cursor()
		# 6. Update the preview 
		self.draw_preview_rectangle()

	
	def on_about(self):
		msg = """ Edit 2D geophysical field.  """
		QMessageBox.about(self, "About", msg.strip())
	
	
	def save_data(self):
		"""
		Saves the data to the netCDF4 file from which input data was read in. First, the program asks
		for a variable name to use, but it remembers this variable name for subsequent saves. 
		"""
		if not self.save_var:
			self.save_var, ok = QInputDialog.getText(self, "Saving...", "Enter variable name:",)
			self.save_var = str(self.save_var)
			if (not ok):
				self.statusBar().showMessage('Save cancelled', 2000)
				self.save_var = None
				return

		ncfile = Dataset(self.dc.fname, "a", format="NETCDF4")
		if not self.save_var in ncfile.variables.keys():
			dvar = ncfile.createVariable(self.save_var, 'f4', (self.dc.lat_var, self.dc.lon_var), zlib=True)
			dvar.units = "km"
		else:
			dvar = ncfile.variables[self.save_var]
		dvar[:,:] = self.dc.data/self.dc.scale
		ncfile.close()
		self.unsaved_changes_exist = False
		self.statusBar().showMessage('Saved to variable: %s' % self.save_var, 2000)

	
	
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
			tip='About TopoEditor')
		self.add_actions(self.help_menu, (about_action,))


	def closeEvent(self, event):
		if self.unsaved_changes_exist:
			reply = QMessageBox.question(self, 'Message', "There are unsaved changes. Are you sure you want to quit?", 
					QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

			if reply == QMessageBox.Yes:
				event.accept()
			else:
				event.ignore()

	
	
def main():
	app = QApplication([])   # Create an application
	app.setWindowIcon(QIcon('Resources/topoicon.png'))

	parser = argparse.ArgumentParser(description='TopoEditor', add_help=False)
	parser.add_argument('fname', nargs=1, type=str, help='name of the netcdf4 data file')
	parser.add_argument('var',   nargs=1, type=str, help='name of the variable in the netcdf4 file')
	parser.add_argument('-s',    nargs=1, type=int, help='size of the view in number of pixels', default=[60])
	parser.add_argument('--scale', nargs=1, type=float, help='multiplicative scaling factor for the data', default=[1.0])
	args = parser.parse_args()

	mw = TopoEditor(args.fname[0], args.var[0], dwx=args.s[0], dwy=args.s[0], scale=args.scale[0])
	mw.show()     # Render the window
	mw.raise_()   # Bring the PyQt4 window to the front
	app.exec_()   # Run the application loop


if __name__ == "__main__":
	main()
