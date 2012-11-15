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
from functions import *
from layer import *

class Facet:
    def __init__(self):
        self.normal = Point()
        self.points = (Point(), Point(), Point())

    def __str__(self):
        s = 'normal: ' + str(self.normal)
        s += ' points:'
        for p in self.points:
            s += str(p)
        return s
    
    def change_direction(self, direction):
        if direction == "+X":
            for p in self.points:
                p.x, p.z = p.z, p.x
        elif direction == "-X":
            for p in self.points:
                p.x, p.z = p.z, -p.x
        elif direction == "+Y":
            for p in self.points:
                p.y, p.z = p.z, p.y
        elif direction == "-Y":
            for p in self.points:
                p.y, p.z = p.z, -p.y
        elif direction == '-Z':
            for p in self.points:
                p.z = -p.z
        elif direction == '+Z':
            pass
        else:
            assert 0

    def intersect(self, z):
        L1 = [True for p in self.points if p.z > z]
        L2 = [True for p in self.points if p.z < z]
        if len(L1) == 3 or len(L2) == 3:
            return (NOT_INTERSECTED, None)
        
        L1 = []
        L2 = []
        for i in range(3):
            p = self.points[i]
            if equal(p.z, z):
                L1.append(i)
            else:
                L2.append(i)
        
        points = self.points
        n = len(L1)
        if n == 0:
            line = self.intersect_0_vertex(points, z)
            code = INTERSECTED
        elif n == 1:
            i1 = L2[0]
            i2 = L2[1]
            p1 = points[i1]
            p2 = points[i2]
            if is_intersected(p1, p2, z):
                line = self.intersect_1_vertex(points[L1[0]], p1, p2, z)
                code = INTERSECTED
            else:
                line = None
                code = NOT_INTERSECTED
        elif n == 2 or n == 3:
            code = REDO
            line = None
        
        return (code, line)

    def intersect_0_vertex(self, points, z):
        L = []
        for i in range(3):
            next = (i + 1) % 3
            p1 = points[i]
            p2 = points[next]
            if is_intersected(p1, p2, z):
                p = calc_intersected_point(p1, p2, z)
                L.append(p)
        
        assert len(L) == 2
        return Line(L[0], L[1])

    def intersect_1_vertex(self, p1, p2, p3, z):
        p = calc_intersected_point(p2, p3, z)
        return Line(p1, p)