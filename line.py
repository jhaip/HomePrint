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
from functions import *
from facet import *
from layer import *

class Line:
    def __init__(self, p1=Point(), p2=Point()):
        self.p1 = p1
        self.p2 = p2

    def __str__(self):
        return str(self.p1) + " -> " + str(self.p2)

    def length(self):
        dx = self.p1.x - self.p2.x
        dy = self.p1.y - self.p2.y
        dz = self.p1.z - self.p2.z
        sum = dx * dx + dy * dy + dz * dz
        return math.sqrt(sum)
    
    def slope(self):
        diffy = self.p2.y - self.p1.y 
        diffx = self.p2.x - self.p1.x
        
        if equal(diffx, 0.0):
            return sys.maxint
        else:
            k = diffy / diffx
            return k