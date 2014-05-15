import os
import glob
import shutil
import re
import time
import math

import pyfits
import numpy as np
import pylab
import matplotlib
import matplotlib.pyplot as plt

USE_PLOT_GUI=False

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.ticker as mticker

import threedhst

def plot_init(square=True, xs=6, aspect=1, left=0.22, bottom=0.11, right=0.02, top=0.02, wspace=0.2, hspace=0.02, fontsize=10, NO_GUI=False, use_tex=False):
    """
    Wrapper for generating a plot window, contains input parameters for setting the 
    full window geometry and also handles toggling the GUI/interactive backend.
    
    NO_GUI should be set to True if your session has no X11 connection.    
    """
    import unicorn
    
    #### If logged in to an external machine ("uni"), don't use GUI plotter
    if unicorn.hostname().startswith('uni') | NO_GUI:
        unicorn.plotting.USE_PLOT_GUI = False
    else:
        unicorn.plotting.USE_PLOT_GUI = True
    
    # plt.rcParams['font.family'] = 'serif'
    # plt.rcParams['font.serif'] = ['Times']
    plt.rcParams['patch.edgecolor'] = 'None'
    plt.rcParams['font.size'] = fontsize

    plt.rcParams['image.origin'] = 'lower'
    plt.rcParams['image.interpolation'] = 'nearest'

    if use_tex:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Times'
        
    if square:
        #xs=5
        lrbt = np.array([left,right,bottom,top])*5./xs     
        ys = (1-lrbt[1]-lrbt[0])/(1-lrbt[3]-lrbt[2])*xs*aspect
        lrbt[[2,3]] /= aspect

        if USE_PLOT_GUI:
            fig = plt.figure(figsize=(xs,ys), dpi=100)
        else:
            fig = Figure(figsize=(xs,ys), dpi=100)
            
        fig.subplots_adjust(left=lrbt[0], bottom=lrbt[2], right=1-lrbt[1], top=1-lrbt[3], wspace=wspace, hspace=hspace)

    else:
        if USE_PLOT_GUI:
            fig = plt.figure(figsize=(7,5), dpi=100)
        else:
            fig = Figure(figsize=(7,5), dpi=100)
            
        fig.subplots_adjust(wspace=wspace, hspace=hspace,left=0.10,
                        bottom=0.10,right=0.99,top=0.97)        
    
    return fig
    
def savefig(fig, filename='figure.png', no_tex=True, dpi=100, increment=False, transparent=False):
    """
    Wrapper around the `savefig` method to handle the two different backends, set whether or not
    an X11/interactive connection is available.
    
    If `increment` is set and the output filename exists, add an integer to
    the filename to avoid overwriting the current figure.
    """
    
    if increment:
        if os.path.exists(filename):
            spl = filename.split('.')
            root, ext = '.'.join(spl[:-1]), spl[-1]
            saved = glob.glob('%s.[0-9]*.%s' %(root, ext))
            filename = '%s.%03d.%s' %(root, len(saved)+1, ext)
            print 'Save %s' %(filename)
            
    if USE_PLOT_GUI:
        fig.savefig(filename,dpi=dpi,transparent=transparent)
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(filename, dpi=dpi, transparent=transparent)
    #
    if no_tex:
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'sans-serif'

class MyLocator(mticker.MaxNLocator):
    """
    Set maximum number of ticks, from
    http://matplotlib.sourceforge.net/examples/pylab_examples/finance_work2.html
    """
    def __init__(self, *args, **kwargs):
        mticker.MaxNLocator.__init__(self, *args, **kwargs)
    
    def __call__(self, *args, **kwargs):
        return mticker.MaxNLocator.__call__(self, *args, **kwargs)


def scatter_annotate(x, y, labels, xtol=None, ytol=None, ax=None,*args, **kwargs):
    if ax is None:
        axi = plt
    else:
        axi = ax
        
    if len(labels) != len(x):
        threedhst.showMessage('`x`, `y`, and `labels` inputs must be same size', warn=True)
        return False
    
    if not isinstance(labels, list):
        labels = list(labels)
        
    plt.scatter(x, y, *args, **kwargs)
    af = AnnoteFinder(x, y, labels, xtol=xtol, ytol=ytol, axis=ax)
    plt.connect('button_press_event', af)
    
    return af
    
    
#  Code taken from http://www.scipy.org/Cookbook/Matplotlib/Interactive_Plotting for annotating
#  a plot label with some sort of label per point.

def linkAnnotationFinders(afs):
  for i in range(len(afs)):
    allButSelfAfs = afs[:i]+afs[i+1:]
    afs[i].links.extend(allButSelfAfs)

class AnnoteFinder:
  """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.
    
  Register this function like this:
    
  scatter(xdata, ydata)
  af = AnnoteFinder(xdata, ydata, annotes)
  connect('button_press_event', af)
  """

  def __init__(self, xdata, ydata, annotes, axis=None, xtol=None, ytol=None):
    self.data = zip(xdata, ydata, annotes)
    if xtol is None:
      # xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
      xtol = (max(xdata) - min(xdata))/40.
    if ytol is None:
      # ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
      ytol = (max(ydata) - min(ydata))/40.
    
    self.xtol = xtol
    self.ytol = ytol
    if axis is None:
      self.axis = pylab.gca()
    else:
      self.axis= axis
    self.drawnAnnotations = {}
    self.links = []

  def distance(self, x1, x2, y1, y2):
    """
    return the distance between two points
    """
    return math.hypot(x1 - x2, y1 - y2)

  def __call__(self, event):
    tb = pylab.get_current_fig_manager().toolbar
    if not (tb.mode == ''):
        return False
        
    if (event.inaxes):
      clickX = event.xdata
      clickY = event.ydata
      if self.axis is None or self.axis==event.inaxes:
        annotes = []
        for x,y,a in self.data:
          if  clickX-self.xtol < x < clickX+self.xtol and  clickY-self.ytol < y < clickY+self.ytol :
            annotes.append((self.distance(x,clickX,y,clickY),x,y, a) )
        if annotes:
          annotes.sort()
          distance, x, y, annote = annotes[0]
          self.drawAnnote(event.inaxes, x, y, annote)
          for l in self.links:
            l.drawSpecificAnnote(annote)

  def drawAnnote(self, axis, x, y, annote):
    """
    Draw the annotation on the plot
    """
    if (x,y) in self.drawnAnnotations:
      markers = self.drawnAnnotations[(x,y)]
      for m in markers:
        m.set_visible(not m.get_visible())
      self.axis.figure.canvas.draw()
    else:
      t = axis.text(x,y, "%s" %(annote), )
      m = axis.scatter([x],[y], marker='s', c='r', zorder=100)
      self.drawnAnnotations[(x,y)] =(t,m)
      self.axis.figure.canvas.draw()

  def drawSpecificAnnote(self, annote):
    annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
    for x,y,a in annotesToDraw:
      self.drawAnnote(self.axis, x, y, a)