import wx
import os
import sys
import string
import copy
import time
import logging
import pprint
import math
import random
import thread
import Queue

try:
    import psyco
    psyco.full()
except ImportError, e:
    pass

try:
    from wx import glcanvas
except ImportError, e:
    print e
    sys.exit()

try:
    from OpenGL.GL import *
    from OpenGL.GLUT import *
except ImportError, e:
    print e
    sys.exit()

ERROR = 2
REDO = 3
LAYER = 4
NOT_LAYER = 5
INTERSECTED = 6
NOT_INTERSECTED = 7
SCANLINE = 8
NOT_SCANLINE = 9
LIMIT = 1e-8

from cadmodel import *
from error import *
from point import *
from line import *
from facet import *
from layer import *

def equal(f1, f2):
    if abs(f1 - f2) < LIMIT:
        return True
    else:
        return False
        
def intersect(x1, y1, x2, y2, x):
    ''' compute y'''
    y = (y2 - y1) / (x2 - x1) * (x - x1) + y1
    return y

def is_intersected(p1, p2, z):
    if (p1.z - z) * (p2.z - z) <= 0.0:
        return True
    else:
        return False

def calc_intersected_point(p1, p2, z):
    x1 = p1.x
    y1 = p1.y
    z1 = p1.z

    x2 = p2.x
    y2 = p2.y
    z2 = p2.z
    
    x = intersect(z1, x1, z2, x2, z)
    y = intersect(z1, y1, z2, y2, z)
    p = Point(x, y, z)
    return p

def writeline(line, f):
    #print >> f, '<line>'
    for p in (line.p1, line.p2):
        print >> f, '<point>',
        print >> f, '<x>', p.x, '</x>',
        print >> f, '<y>', p.y, '</y>',
        print >> f, '<z>', p.z, '</z>',
        print >> f, '</point>'
    #print >> f, '</line>'        