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
from facet import *

class Layer:
    colors = ([1, 0, 1], [0, 1, 1], [1, 1, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, 1])

    def __init__(self, z, pitch):
        self.lines = []
        self.z = z
        self.pitch = pitch

    def empty(self):
        return len(self.lines) == 0

    def create_gllist(self):
        self.layerListId = 1001
        glNewList(self.layerListId, GL_COMPILE)
        
        glBegin(GL_LINES)
        for chunk in self.chunks:
            r = random.random()
            g = random.random()
            b = random.random()
            glColor(r, g, b)
 
            for line in chunk:
                for p in [line.p1, line.p2]:
                    glVertex3f(p.x, p.y, p.z)
        
        glColor(1, 1, 0)
        for loop in self.loops:
            for line in loop:
                for p in [line.p1, line.p2]:
                    glVertex3f(p.x, p.y, p.z)
        
        glEnd()

        glBegin(GL_POINTS)
        glColor(1, 0, 1)
        for loop in self.loops:
            counter = 0
            for line in loop:
                #for p in [line.p1, line.p2]:
                if counter == 0:
                    glColor(0,1,1)
                elif counter == 1:
                    glColor(1,0,0)
                else:
                    glColor(1,0,1)
                glVertex3f(line.p1.x, line.p1.y, line.p1.z)
                counter += 1
        glEnd()

        glEndList()
        return self.layerListId

    def set_lines(self, lines):
        self.lines = lines
        ok = self.createLoops()
        if not ok:
            return False
        
        self.calc_dimension()             
        #self.create_scanlines()
        #self.create_chunks()
        self.scanlines = [] ##
        self.chunks = [] ##

        return True

    def createLoops(self):
        lines = self.lines

        self.loops = []
        while len(lines) != 0:
            loop = []
            line = lines.pop()
            loop.append(line)
            
            start = line.p1
            p2 = line.p2
            while True:
                found = False
                for aline in lines:
                    if p2 == aline.p1:
                        p1 = aline.p1
                        p2 = aline.p2
                        found = True
                        break
                    elif p2 == aline.p2:
                        p1 = aline.p2
                        p2 = aline.p1
                        found = True
                        break

                if found:        
                    lines.remove(aline)
                    loop.append(Line(p1, p2))
                    if p2 == start:
                        break
                else:
                    print 'error: loop is not found'
                    return False
            
            self.move_lines(loop)
            nloop = self.merge_lines(loop)
            self.loops.append(nloop)
        
        return True                
    
    def move_lines(self, loop):
        tail = loop[-1]
        k1 = tail.slope()
        head = loop[0]
        k2 = head.slope()
        rm_list = []
        if equal(k1, k2):
            for aline in loop:
                k = aline.slope()
                if equal(k, k1):
                    rm_list.append(aline)
                else:
                    break
            
            for it in rm_list:
                loop.remove(it)
            
            loop.extend(rm_list)
        
        k1 = loop[0].slope()
        k2 = loop[-1].slope()
        assert not equal(k1, k2)

    def merge_lines(self, loop):
        nloop = []
        while len(loop) != 0:
            line = loop.pop(0) 
            k1 = line.slope()
            p1 = line.p1
            p2 = line.p2
            rm_list = []            
            for aline in loop:
                k2 = aline.slope()
                if equal(k1, k2):
                    p2 = aline.p2
                    rm_list.append(aline)
                else:
                    p2 = aline.p1
                    break
            
            for it in rm_list:
                loop.remove(it)
            nloop.append(Line(p1, p2))
        
        return nloop

    def calc_dimension(self):
        ylist = []
        for loop in self.loops:
            for line in loop:
                ylist.append(line.p1.y)
                ylist.append(line.p2.y)
        self.miny = min(ylist)                
        self.maxy = max(ylist)

    def intersect(self, y, line, loop):
        y1 = line.p1.y
        y2 = line.p2.y
        if self.is_intersected(y1, y2, y):
            count = 0
            if equal(y, y1):
                count += 1
                p = line.p1

            if equal(y, y2):
                count += 1
                p = line.p2
            
            if count == 0:
                x = self.intersect_0(y, line)
                code = INTERSECTED
            elif count == 1:
                if self.is_peak(y, p, line, loop):
                    code = NOT_INTERSECTED
                    x = None
                else:
                    code = INTERSECTED
                    x = p.x
            elif count == 2:
                code = REDO
                x = None
        else:
            code = NOT_INTERSECTED 
            x = None
        
        return (code, x)

    def intersect_0(self, y, line):
        x1 = line.p1.x
        y1 = line.p1.y
        x2 = line.p2.x
        y2 = line.p2.y
        
        if equal(x1, x2):
            x = x1
            return x
        else:
           x = (y -  y1) * (x2 - x1) / (y2 - y1) + x1
           return x
    
    def is_peak(self, y, point, line, loop):
        L = []
        for it in loop:
            if point == it.p1:
                L.append(it.p2)
            elif point == it.p2:
                L.append(it.p1)
        
        val = (L[0].y - y) * (L[1].y - y)
        if val > 0.0:
            return True
        else:
            return False

    def is_intersected(self, y1, y2, y):
        if (y1 - y) * (y2 - y) <= 0.0:
            return True
        else:
            return False
    
    def get_overlap_line(self, line, scanline):
        y2 = scanline[0].p1.y
        y1 = line.p1.y
        
        # Are they adjacent lines?
        distance = abs(y2 - y1)
        if equal(distance, self.pitch) or distance < self.pitch:
            for aline in scanline:
                if aline.p1.x >= line.p2.x or aline.p2.x <= line.p1.x:
                    pass
                else:
                    return aline
        else:
            return False  

    def distance(xi,xii,yi,yii):
        sq1 = (xi-xii)*(xi-xii)
        sq2 = (yi-yii)*(yi-yii)
        return math.sqrt(sq1 + sq2)

    def write(self, f, global_start):
        #print >> f, '<layer id="', self.id, '">'
        len = self.writeloop(f, global_start)
        #self.writechunks(f)
        #print >> f, '</layer>'
        print "LINE LENGTH"
        print len
        print >> f, 'WAIT SEC 30.0'
    
    def writeloop(self, f, global_start):
        #print >> f, '<loops num="', len(self.loops), '">'
        count = 1
        looplength = 0
        for loop in self.loops:
            #print >> f, '<loop id="', count, '">'
            print >> f, 'OUT[4] = TRUE'
            lastline = None
            for line in loop:
                looplength += self.distance(line.p1.x,line.p2.x,line.p1.y,line.p2.y)
                if lastline != None:
                    if (int(line.p1.x) != int(lastline.p2.x)) or (int(line.p1.y) != int(lastline.p2.y)) or (int(line.p1.z) != int(lastline.p2.z)):
                        print >> f, 'LIN {X '+str(int(global_start.x+line.p1.x*1000))+', Y '+str(int(global_start.y+line.p1.y*1000))+', Z '+str(int(global_start.z+line.p1.z*1000))+'} c_vel'
                    print >> f, 'LIN {X '+str(int(global_start.x+line.p2.x*1000))+', Y '+str(int(global_start.y+line.p2.y*1000))+', Z '+str(int(global_start.z+line.p2.z*1000))+'} c_vel'
                else:
                    print >> f, 'LIN {X '+str(int(global_start.x+line.p1.x*1000))+', Y '+str(int(global_start.y+line.p1.y*1000))+', Z '+str(int(global_start.z+line.p1.z*1000))+'} c_vel'
                    print >> f, 'LIN {X '+str(int(global_start.x+line.p2.x*1000))+', Y '+str(int(global_start.y+line.p2.y*1000))+', Z '+str(int(global_start.z+line.p2.z*1000))+'} c_vel'
                #writeline(line, f)
                lastline = line
            count += 1
            print >> f, 'OUT[4] = FALSE'
        #print >> f, '</loops>'
        return looplength

    def writechunks(self, f):
        print >> f, '<chunks num="', len(self.chunks), '">'
        count = 1
        for chunk in self.chunks:
            print >> f, '<chunk id="', count, '">'
            for line in chunk:
                writeline(line, f)
            print >> f, '</chunk>'
            count += 1
        print >> f, '</chunks>'