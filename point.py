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
from line import *
from functions import *
from facet import *
from layer import *

class Point:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        s = '(%f, %f, %f) ' % (self.x, self.y, self.z)
        return s

    def __eq__(self, other):
        return equal(self.x, other.x) and equal(self.y, other.y) and equal(self.z, other.z)

    def __cmp__(self, other):
        if self == other:
            return 0
        elif self.x < other.x or self.y < other.y or self.z < other.z:
            return -1
        else:
            return 1
    
    def __hash__(self):
        s = '%.6f %.6f %.6f' % (self.x, self.y, self.z)
        return hash(s)