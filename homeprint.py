#!/usr/bin/env python 
#-----------------------------------------------------------------------------
# Author     : Zhigang Liu
# Date       : Jan 2009
# Email      : zgliu71@gmail.com
# License    : General Public License 2 (GPL2) 
# Description: Slice STL CAD file layer by layer
#-----------------------------------------------------------------------------

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

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

def equal(f1, f2):
    if abs(f1 - f2) < LIMIT:
        return True
    else:
        return False

class EndFileException(Exception):
    def __init__(self, args=None):
        self.args = args

class FormatError(Exception):
    def __init__(self, value=None):
        self.value = value
    
    def __str__(self):
        return 'FormatError:' + self.value

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

def writeline(line, f):
    #print >> f, '<line>'
    for p in (line.p1, line.p2):
        print >> f, '<point>',
        print >> f, '<x>', p.x, '</x>',
        print >> f, '<y>', p.y, '</y>',
        print >> f, '<z>', p.z, '</z>',
        print >> f, '</point>'
    #print >> f, '</line>'        

class CadModel:
    def __init__(self):
        self.init_logger()
        self.loaded = False
        self.curr_layer = -1
        self.sliced = False
        self.dimension = {}
        self.old_dimension = {}
        self.scale = 1
    
    def next_layer(self):
        n = len(self.layers)
        self.curr_layer = (self.curr_layer + 1) % len(self.layers)
    
    def prev_layer(self):
        n = len(self.layers)
        self.curr_layer -= 1
        if self.curr_layer == -1:
            self.curr_layer = len(self.layers) -1

    def get_curr_layer(self):
        return self.layers[self.curr_layer]

    def init_logger(self):
        #self.logger = logging.getLogger(self.__class__.__name__)
        self.logger = logging.getLogger("cadmodel")
        self.logger.setLevel(logging.DEBUG)
        h = logging.StreamHandler()
        h.setLevel(logging.DEBUG)
        f = logging.Formatter("%(levelname)s %(filename)s:%(lineno)d %(message)s")
        h.setFormatter(f)
        self.logger.addHandler(h)
    
    def get_line(self, f):
        line = f.readline()
        if not line:
            raise EndFileException, 'end of file'
        return line.strip()

    def get_normal(self, f):
        line = self.get_line(f)
        items = line.split()
        num_items = len(items)
        if num_items != 5:
            if num_items >=1 and items[0] == "endsolid":
                self.loaded = True
                raise EndFileException, 'endfile'
            else:
                self.logger.error(line)
                raise FormatError, line
        
        if items[0] != 'facet' and items[1] != 'normal':
            self.logger.error(line)
            raise FormatError, line
        
        L = map(lambda x: float(x), items[2:])
        normal = Point(L[0], L[1], L[2])
        return normal

    def get_outer_loop(self, f):
        line = self.get_line(f)
        if line != "outer loop":
            self.logger.error(line)
            raise FormatError, line

    def get_vertex(self, f):
        points = []
        for i in range(3):
            line = self.get_line(f)
            items = line.split()
            no = len(items)
            if no != 4:
                self.logger.error(line)
                raise FormatError, line
            if items[0] != 'vertex':
                self.logger.error(line)
                raise FormatError, line

            L = map(lambda x: float(x), items[1:])
            point = Point(L[0], L[1], L[2])
            points.append(point)
        return points
    
    def get_end_loop(self, f):
        line = self.get_line(f) 
        if line != 'endloop':
            self.logger.error(line)
            raise FormatError, line
    
    def get_end_facet(self, f):
        line = self.get_line(f)
        if line != 'endfacet':
            self.logger.error(line)
            raise FormatError, line

    def get_facet(self, f):
        normal = self.get_normal(f)   
        self.get_outer_loop(f)
        points = self.get_vertex(f)
        facet = Facet()
        facet.normal = normal
        facet.points = points
        self.get_end_loop(f)
        self.get_end_facet(f)
        return facet
    
    def get_solid_line(self, f):
        ''' Read the first line'''
        line = self.get_line(f)
        items = line.split()
        no = len(items)
        if no >= 2 and items[0] == 'solid':
            self.modelName = items[1]
        else:
            self.logger.error(line)
            raise FormatError, line
    
    def calc_dimension(self):
        if self.loaded:
            xlist = []
            ylist = []
            zlist = []
            for facet in self.facets:
                for p in facet.points:
                    xlist.append(p.x)
                    ylist.append(p.y)
                    zlist.append(p.z)
            self.minx = min(xlist)
            self.maxx = max(xlist)
            self.miny = min(ylist)
            self.maxy = max(ylist)
            self.minz = min(zlist)
            self.maxz = max(zlist)
            
            self.xsize = self.maxx - self.minx
            self.ysize = self.maxy - self.miny
            self.zsize = self.maxz - self.minz

            self.diameter = math.sqrt(self.xsize * self.xsize + self.ysize * self.ysize + self.zsize * self.zsize)

            # Center
            self.xcenter = (self.minx + self.maxx) / 2
            self.ycenter = (self.miny + self.maxy) / 2
            self.zcenter = (self.minz + self.maxz) / 2

    def open(self, filename):
        start = time.time()
        try:
            f = open(filename) 
        except IOError, e:
            print e
            return False
        
        try:
            self.get_solid_line(f)
            self.facets = [] 
            while True:
                facet = self.get_facet(f)
                self.facets.append(facet)
        except EndFileException, e:
            pass
        except FormatError, e:
            print e
            return False
        
        if self.loaded:
            self.calc_dimension()
            self.logger.debug("no of facets:" + str(len(self.facets)))
            self.oldfacets = copy.deepcopy(self.facets)
            self.sliced = False
            self.set_old_dimension()
            cpu = '%.1f' % (time.time() - start)
            
            print 'open cpu', cpu, 'secs'
            return True
        else:
            return False

    def save(self, filename):
        f = open(filename, 'w')

        print >> f, '&ACCESS RVP'
        print >> f, '&REL 25'
        print >> f, 'DEF '+filename+'()'
        print >> f, ''
        print >> f, ';FOLD INI'
        print >> f, '  ;FOLD BASISTECH INI'
        print >> f, '    GLOBAL INTERRUPT DECL 3 WHEN $STOPMESS==TRUE DO IR_STOPM ( )'
        print >> f, '    INTERRUPT ON 3'
        print >> f, '    BAS (#INITMOV,0 )'
        print >> f, '  ;ENDFOLD (BASISTECH INI)'
        print >> f, '  ;FOLD USER INI'
        print >> f, '    ;Make your modifications here'
        print >> f, ''
        print >> f, '  ;ENDFOLD (USER INI)'
        print >> f, ';ENDFOLD (INI)'
        print >> f, ''
        print >> f, ';FOLD PTP HOME  Vel= 100 % DEFAULT;%{PE}%R 5.4.37,%MKUKATPBASIS,%CMOVE,%VPTP,%P 1:PTP, 2:HOME, 3:, 5:100, 7:DEFAULT'
        print >> f, '$BWDSTART=FALSE'
        print >> f, 'PDAT_ACT=PDEFAULT'
        print >> f, 'FDAT_ACT=FHOME'
        print >> f, 'BAS(#PTP_PARAMS,100)'
        print >> f, '$H_POS=XHOME'
        print >> f, 'PTP XHOME'
        print >> f, ';ENDFOLD'
        print >> f, ''
        print >> f, '$VEL_AXIS[1]=20'
        print >> f, '$VEL_AXIS[2]=20'
        print >> f, '$VEL_AXIS[3]=20'
        print >> f, '$VEL_AXIS[4]=20'
        print >> f, '$VEL_AXIS[5]=20'
        print >> f, '$VEL_AXIS[6]=20'

        print >> f, ''
        print >> f, '$vel.cp= '+str(self.speed)
        print >> f, '$apo.cvel= 95'
        print >> f, ''

        for layer in self.layers:
            layer.write(f, self.global_start)

        print >> f, 'END'

    def slice(self, para):
        self.sliced = False
        self.height = float(para["height"])/1000.0
        self.pitch = float(para["pitch"])
        self.speed = float(para["speed"])
        #self.fast = float(para["fast"])
        self.direction = para["direction"]
        #self.scale = float(para["scale"])
        self.global_start = Point(float(para["global_start_x"]),float(para["global_start_y"]),float(para["global_start_z"]))
        
        #self.scale_model(self.scale)
        self.change_direction(self.direction)
        self.calc_dimension()
        self.create_layers()
        self.set_new_dimension()
        self.change_direction('+Z')
        if len(self.layers) > 0:
            self.sliced = True
            self.curr_layer = 0
            return True
        else:
            self.sliced = False
            return False
    
    def set_old_dimension(self):
        self.old_dimension["x"] = self.xsize
        self.old_dimension["y"] = self.ysize
        self.old_dimension["z"] = self.zsize
        self.dimension["x"] = str(self.xsize)
        self.dimension["y"] = str(self.ysize)
        self.dimension["z"] = str(self.zsize)
        self.dimension["factor"] = str(self.scale)

    def set_new_dimension(self):
        self.dimension["x"] = str(self.xsize)
        self.dimension["y"] = str(self.ysize)
        self.dimension["z"] = str(self.zsize)
        self.dimension["factor"] = str(self.scale)

    def scale_model(self, factor):
        self.facets = []
        for facet in self.oldfacets:
            nfacet = copy.deepcopy(facet)
            for p in nfacet.points:
                p.x *= factor
                p.y *= factor
                p.z *= factor
            self.facets.append(nfacet)
    
    def change_direction(self, direction):
        for facet in self.facets:
            print "changing direction"
            facet.change_direction(direction)
    
    def create_layers(self):
        start = time.time()
        self.layers = []
        z = self.minz + self.height
        lastz = self.minz
        count = 0

        no = (self.maxz - self.minz) / self.height
        no = int(no)
        self.queue.put(no)
        while z >= self.minz and z <= self.maxz:
            print 'HIT'
            code, layer = self.create_one_layer(z)
            
            if code == LAYER:
                count += 1
                layer.id = count
                self.layers.append(layer)
                
                lastz = z
                z += self.height
                self.queue.put(count)
                print 'layer', count, '/', no
            elif code == ERROR:
                print 'layer error'
                break
            elif code == REDO:
                z = z - self.height * 0.01
                if z < lastz:
                    break
                print 'recreate layer'
            elif code == NOT_LAYER:
                lastz = z
                z += self.height
                print 'not layer'
           
        self.queue.put("done")                
        print 'no of layers:', len(self.layers)                
        cpu = '%.1f' % (time.time() - start)
        print 'slice cpu', cpu,'secs'
    
    def create_one_layer(self, z):
        layer = Layer(z, self.pitch)
        lines = []
        for facet in self.facets:
            print 'facet check'
            code, line = facet.intersect(z) 
            if code == REDO:
                return (REDO, None)
            elif code == INTERSECTED:
                lines.append(line)
        
        print 'cleanup stuff'
        if len(lines) != 0:
            ok = layer.set_lines(lines)
            if ok:
                return (LAYER, layer)
            else:
                return (ERROR, None)
        else:
            return (NOT_LAYER, None)
    
    def create_gl_model_list(self):
        wireframe = True
        self.model_list_id = 1000
        glNewList(self.model_list_id, GL_COMPILE)
        if self.loaded:
            for facet in self.facets:
                if wireframe == False:
                    normal = facet.normal
                    glNormal3f(normal.x, normal.y, normal.z)
                    glColor(0.5,0.5,0.5)
                    glBegin(GL_TRIANGLES)
                else:
                    glColor(1,1,1)
                    glBegin(GL_LINE_LOOP)
                for p in facet.points:
                    glVertex3f(p.x, p.y, p.z)
                glEnd()
        glEndList()

    def create_gl_layer_list(self):
        assert self.sliced
        layer = self.get_curr_layer()
        return layer.create_gllist()

class PathCanvas(glcanvas.GLCanvas):
    def __init__(self, parent, cadmodel):
        glcanvas.GLCanvas.__init__(self, parent, -1)

        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.cadmodel = cadmodel

    def OnEraseBackground(self, event):
        pass

    def OnSize(self, event):
        if self.GetContext():
            self.SetCurrent()
            size = self.GetClientSize()
            glViewport(0, 0, size.width, size.height)
        self.Refresh()
        event.Skip()

    def OnPaint(self, event):
        dc = wx.PaintDC(self)
        self.SetCurrent()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.show_path()
        self.SwapBuffers()

    def setup_projection(self):
        diameter = self.cadmodel.diameter
        size = self.GetClientSize()
        w = size.width
        h = size.height
        
        half = diameter / 2
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()

        if w <= h:
            factor = float(h) / w
            left = -half
            right = half
            bottom = -half * factor
            top = half * factor
        else:
            factor = float(w) / h
            left  = -half * factor 
            right = half * factor
            bottom = -half
            top = half
        near = 0
        far = diameter * 2
        glOrtho(left, right, bottom, top, near, far)

        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE)
        glPointSize(5.0)

    def show_path(self):
        if self.cadmodel.sliced:
            self.setup_projection()
            glMatrixMode(GL_MODELVIEW)
            glLoadIdentity()
            layer = self.cadmodel.get_curr_layer()
            z = layer.z
            glTranslatef(-self.cadmodel.xcenter, -self.cadmodel.ycenter, -z)
            layer_id = self.cadmodel.create_gl_layer_list()
            glCallList(layer_id)
            
class ModelCanvas(glcanvas.GLCanvas):

    def __init__(self, parent, cadmodel):
        glcanvas.GLCanvas.__init__(self, parent, -1)
        self.init = False
        self.cadmodel = cadmodel
        self.lastx = self.x = 30
        self.lasty = self.y = 30
        self.xangle = 0
        self.yangle = 0
        self.zoom = 1

        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnMouseDown)
        self.Bind(wx.EVT_LEFT_UP, self.OnMouseUp)
        self.Bind(wx.EVT_MOTION, self.OnMouseMotion)
        self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)

    def OnEraseBackground(self, event):
        pass # Do nothing, to avoid flashing on MSW.

    def OnPaint(self, event):
        dc = wx.PaintDC(self)
        self.SetCurrent()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.show_model()
        #self.show_path()
        self.show_grid()
        self.show_axis()
        self.SwapBuffers()

    def show_grid(self):
        glLineWidth(1.)
        glBegin(GL_LINES)
        cellsize = 0.1
        gridw = 10

        gridh = 10
        x = -gridw*cellsize
        while x <= gridw*cellsize:
            y = -gridh*cellsize
            while y <= gridh*cellsize:
                #glColor(0,0,1)
                glColor(0,0.5,0.5)
                if x == 0:
                    glColor(1,1,1)
                glVertex3f(x,-gridh*cellsize,0)
                glVertex3f(x,gridh*cellsize,0)
                #glColor(1,0,1)
                if y == 0:
                    glColor(1,1,1)
                glVertex3f(-gridw*cellsize,y,0)
                glVertex3f(gridw*cellsize,y,0)

                y += cellsize
            x += cellsize

        glEnd()
    
    #def show_path(self):
    #    if self.cadmodel.sliced:
    #        layer_id = self.cadmodel.create_gl_layer_list()
    #        glCallList(layer_id)

    def show_model(self):
        if not self.cadmodel.loaded:
            return
        
        #self.setup_gl_context()
        self.setup_projection()
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

        glRotate(45,0,0,1)
         
        glTranslatef(0, 0, -self.cadmodel.diameter)
        # Rotate model
        glRotatef(self.xangle, 1, 0, 0)
        glRotatef(self.yangle, 0, 1, 0)
        
        # Move model to origin
        #glTranslatef(-self.cadmodel.xcenter, -self.cadmodel.ycenter, -self.cadmodel.zcenter)
        
        # Move model to origin, at z = 0
        #glTranslatef(-self.cadmodel.xcenter, -self.cadmodel.ycenter, 0)

        glCallList(self.cadmodel.model_list_id)

        if self.cadmodel.sliced:
            layer_id = self.cadmodel.create_gl_layer_list()
            glCallList(layer_id)

        # Model coordinate system back to (0,0,0)
        #glTranslatef(self.cadmodel.xcenter, self.cadmodel.ycenter, 0)

    def show_axis(self):
        glLineWidth(3.)
        glBegin(GL_LINES)
        glColor(1,0,0)
        glVertex3f(0.,0.,0.)
        glVertex3f(3.,0.,0.)
        glColor(0,1,0)
        glVertex3f(0.,0.,0.)
        glVertex3f(0.,3.,0.)
        glColor(0,0,1)
        glVertex3f(0.,0.,0.)
        glVertex3f(0.,0.,3.)
        glEnd()
        glLineWidth(2.)
        glTranslatef(0.,0.,4.)
        glutStrokeCharacter(GLUT_STROKE_ROMAN, ord('Z'));
        glTranslatef(0.,0.,-4.)

    def OnMouseDown(self, evt):
        self.CaptureMouse()
        self.x, self.y = self.lastx, self.lasty = evt.GetPosition()

    def OnMouseUp(self, evt):
        if self.HasCapture():
            self.ReleaseMouse()

    def OnMouseMotion(self, evt):
        if evt.Dragging() and evt.LeftIsDown():
            self.lastx, self.lasty = self.x, self.y
            self.x, self.y = evt.GetPosition()

            self.xangle += (self.y - self.lasty)
            self.yangle += (self.x - self.lastx)
            self.Refresh(False)

    def OnMouseWheel(self, evt):
        self.zoom += evt.GetWheelRotation()*0.01
        print self.zoom

    def create_model(self):
        self.xangle = 0
        self.yangle = 0
        self.SetCurrent()
        
        if not self.init:
            self.setup_gl_context()
            self.init =  True
        self.cadmodel.create_gl_model_list()
        self.Refresh()

    def OnSize(self, event):
        if self.GetContext():
            self.SetCurrent()
            self.setup_viewport()
        self.Refresh()
        event.Skip()
    
    def setup_viewport(self):
        size = self.GetClientSize()
        glViewport(0, 0, size.width, size.height)

    def setup_projection(self):
        maxlen = self.cadmodel.diameter
        size = self.GetClientSize()
        w = size.width
        h = size.height
        
        half = maxlen / 2
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()

        if w <= h:
            factor = float(h) / w
            left = -half
            right = half
            bottom = -half * factor
            top = half * factor
        else:
            factor = float(w) / h
            left  = -half * factor 
            right = half * factor
            bottom = -half
            top = half
        near = -maxlen * 4
        far = maxlen * 4
        zoom = 2.
        glOrtho(left*zoom, right*zoom, bottom*zoom, top*zoom, near*zoom, far*zoom)    

    def setup_gl_context(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)

        ambient_light = [0,0,0,1]#[0.2, 0.2, 0.2, 1.0]
        diffuse_light = [1,1,1,1]#[0.8, 0.8, 0.8, 1.0]
        specular_light = [1,1,1,1]#[0.5, 0.5, 0.5, 1.0]
        position = [-1.5, 1.0, -4.0, 1.0]
        position = [1,1,0]#[-15.0, 30.0, -40.0, 1.0]

        glLightfv(GL_LIGHT0, GL_AMBIENT, ambient_light)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_light)
        glLightfv(GL_LIGHT0, GL_SPECULAR, specular_light)
        glLightfv(GL_LIGHT0, GL_POSITION, position)
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.2, 0.2, 0.2, 1.0])

        mcolor = [ 0.0, 0.0, 0.4, 0.0]
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mcolor)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_CULL_FACE)
        glPolygonMode(GL_BACK, GL_LINE)
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE)
        glEnable(GL_COLOR_MATERIAL)
        glMaterial(GL_FRONT, GL_SHININESS, 50)#96)

        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE)
        glPointSize(7.0)

        glShadeModel(GL_SMOOTH)
        glEnable(GL_POINT_SMOOTH)   

class DimensionPanel(wx.Panel):
    def __init__(self, parent, cadmodel):
        wx.Panel.__init__(self, parent)
        self.cadmodel = cadmodel
        self.ids = {"X":5300,"Y":5302,"Z":5303,"Factor":5304}
        self.handlers = {"X":self.OnDimXChange,"Y":self.OnDimYChange,"Z":self.OnDimZChange,"Factor":self.OnFactorChange}
        self.txt_fields = {}
        self.create_controls()

    def create_controls(self):
        box = wx.StaticBox(self, label="Dimension") 
        sizer = wx.StaticBoxSizer(box, wx.HORIZONTAL)
        self.SetSizer(sizer)
        
        label = "Original"
        items = [("X", "x"), ("Y", "y"), ("Z", "z"),("Factor","factor")]
        s1 = self.create_dimension(label, items)
        sizer.Add(s1, 1, wx.EXPAND|wx.ALL, 2)

        #label = "Scaled"
        #items = [("X", 'newx'), ('Y', 'newy'), ('Z', 'newz')]
        #s2 = self.create_dimension(label, items)
        #sizer.Add(s2, 1, wx.EXPAND|wx.ALL, 2)

    def create_dimension(self, label, items):
        sizer = wx.BoxSizer(wx.VERTICAL) 
        caption = wx.StaticText(self, label=label)
        sizer.Add(caption, 0, wx.ALIGN_CENTER)

        flex = wx.FlexGridSizer(rows=len(items), cols=2, hgap=2, vgap=2)
        for label, key in items:
            lbl_ctrl = wx.StaticText(self, label=label)
            txt_ctrl = wx.TextCtrl(self, id=self.ids[label], value="", size=(70, -1), style=wx.TE_PROCESS_ENTER)
            flex.Add(lbl_ctrl)
            flex.Add(txt_ctrl, 0, wx.EXPAND)
            self.txt_fields[key] = txt_ctrl

        for i in self.ids:
            self.Bind(wx.EVT_TEXT_ENTER, self.handlers[i], id=self.ids[i])

        sizer.Add(flex, 0, wx.EXPAND)
        flex.AddGrowableCol(1, 1)
        return sizer

    def set_values(self, dimension):
        for key in dimension:
            self.txt_fields[key].SetValue(dimension[key])

    def update_dimension(self,factor):
        print "Updating dimensions"
        self.cadmodel.scale = factor
        self.cadmodel.scale_model(factor)
        self.cadmodel.calc_dimension()
        self.cadmodel.set_new_dimension()
        self.set_values(self.cadmodel.dimension)

    def OnDimXChange(self, event):
        v = float(self.txt_fields["x"].GetValue())
        self.update_dimension(v/self.cadmodel.old_dimension["x"])
        event.Skip()

    def OnDimYChange(self, event):
        v = float(self.txt_fields["y"].GetValue())
        self.update_dimension(v/self.cadmodel.old_dimension["y"])
        event.Skip()

    def OnDimZChange(self, event):
        v = float(self.txt_fields["z"].GetValue())
        self.update_dimension(v/self.cadmodel.old_dimension["z"])
        event.Skip()

    def OnFactorChange(self, event):
        f = float(self.txt_fields["factor"].GetValue())
        self.update_dimension(f)
        event.Skip()

class ControlPanel(wx.Panel):
    def __init__(self, parent, cadmodel):
        wx.Panel.__init__(self, parent, -1)
        self.cadmodel = cadmodel
        self.create_controls()

    def create_controls(self):
        mainsizer = wx.BoxSizer(wx.VERTICAL)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        mainsizer.Add(sizer, 1, wx.ALL|wx.EXPAND, 10)
        self.SetSizer(mainsizer)
        
        # Dimension panel
        self.dimensionPanel = DimensionPanel(self, self.cadmodel)
        sizer.Add(self.dimensionPanel, 0, wx.EXPAND|wx.ALIGN_CENTER)
        
        # Slice info panel
        sizer.Add((10,10)) 
        sliceSizer = self.create_slice_info()
        sizer.Add(sliceSizer, 0, wx.EXPAND)

        # Slice Button
        box = wx.BoxSizer(wx.HORIZONTAL)
        box.Add(wx.Button(self, -1, 'Slice'), 1 )
        sizer.Add(box, 0, wx.EXPAND)
        
        # Change Axis
        #box2 = wx.BoxSizer(wx.HORIZONTAL)
        #axis = ["+X", "-X", "+Y", "-Y", "+Z", "-Z"]
        #box2.Add(wx.ComboBox(self, value='+Z', pos=wx.DefaultPosition, size=wx.DefaultSize, choices=axis, style=wx.CB_READONLY), 1)
        #sizer.Add(box2, 0, wx.EXPAND)
        #self.Bind(wx.EVT_COMBOBOX, self.OnSelect)

        # image
        # sizer.AddStretchSpacer()
        # img = wx.Image('cat.jpg', wx.BITMAP_TYPE_ANY)
        # w = img.GetWidth()
        # h = img.GetHeight()
        # factor = 0.6
        # img = img.Scale(w * factor, h * factor)
        # img = img.ConvertToGreyscale()
        # sb = wx.StaticBitmap(self, -1, wx.BitmapFromImage(img), style=wx.RAISED_BORDER)
        # sizer.Add(sb, 0, wx.ALIGN_CENTER_HORIZONTAL)

    #def OnSelect(self, event):
    #    item = event.GetSelection()
    #    axis = ['+X', '-X', '+Y', '-Y', '+Z', '-Z']
    #    print "ON SELECT EVENT"
    #    print axis[item]
    #    print "THAT's WHAT HAPPENED"
    #    self.cadmodel.direction = axis[item]
    #    self.cadmodel.change_direction(axis[item])
    #    self.cadmodel.calc_dimension()
    #    self.cadmodel.set_new_dimension()

    def create_slice_info(self):
        self.txt_fields = {}
        box = wx.StaticBox(self, -1, "Slice Info")
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)

        #items = [("Layer hight", "height"), ("Pitch", "pitch"), ("Speed", "speed"), 
        #         ("Direction", "direction"), ("Num Layers", "nolayer"),
        #         ("Current Layer", "currlayer")]
        items = [("Layer hight", "height"), ("Speed", "speed"), 
                 ("Direction", "direction"), ("Num Layers", "nolayer"),
                 ("Current Layer", "currlayer")]
        flex = wx.FlexGridSizer(rows=len(items), cols=2, hgap=2, vgap=2)
        for label, key in items:
            lbl_ctrl = wx.StaticText(self, label=label)
            txt_ctrl = wx.TextCtrl(self, size=(70, -1), style=wx.TE_READONLY)
            flex.Add(lbl_ctrl)
            flex.Add(txt_ctrl, 0, wx.EXPAND)
            self.txt_fields[key] = txt_ctrl
        flex.AddGrowableCol(1, 1)
        sizer.Add(flex, 1, wx.EXPAND|wx.ALL, 2)
        return sizer

    def set_dimension(self, dimension): 
        self.dimensionPanel.set_values(dimension)

    def set_slice_info(self, info):
        for key in self.txt_fields.keys():
            txt = self.txt_fields[key]
            value = info.get(key, "")
            txt.SetValue(value)
    
    def set_num_layer(self, num_layers):
        self.txt_fields["nolayer"].SetValue(str(num_layers))

    def set_curr_layer(self, curr_layer):
        self.txt_fields["currlayer"].SetValue(str(curr_layer))

class BlackcatFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "HomePrint - MIT Mediated Matter", size=(800, 600))
        #self.slice_parameter = {"height":"1.0", "pitch":"1.0", "speed":"10", "fast":"20", "direction":"+Z", "scale":"1"}
        self.slice_parameter = {"height":"1.0", "speed":"0.3", "pitch": "1.0", "direction":"+Z", "global_start_x":"0", "global_start_y":"0", "global_start_z":"0"}
        self.create_menubar()
        self.create_toolbar()
        self.cadmodel = CadModel()
        self.statusbar = self.CreateStatusBar()
        self.create_panel()
        self.Centre()

    def create_toolbar(self):
        self.ID_SLICE = 1001
        self.ID_NEXT = 2000
        self.ID_PREV = 2001

        toolbar = self.CreateToolBar()
        img_open = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN)
        img_save = wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE)
        img_slice = wx.ArtProvider.GetBitmap(wx.ART_CDROM)
        img_next = wx.ArtProvider.GetBitmap(wx.ART_GO_DOWN)
        img_prev = wx.ArtProvider.GetBitmap(wx.ART_GO_UP)
        img_help = wx.ArtProvider.GetBitmap(wx.ART_HELP, client=wx.ART_TOOLBAR)
        img_quit = wx.ArtProvider.GetBitmap(wx.ART_QUIT)

        toolbar.AddLabelTool(wx.ID_OPEN, 'open', img_open, shortHelp='open file', longHelp='open CAD model')
        toolbar.AddLabelTool(self.ID_SLICE, 'slice', img_slice, shortHelp='slice modal')
        toolbar.AddLabelTool(wx.ID_SAVE, 'save', img_save, shortHelp='save slice info', longHelp='save slice result')
        toolbar.AddLabelTool(self.ID_NEXT, 'next', img_next, shortHelp='next layer')
        toolbar.AddLabelTool(self.ID_PREV, 'prev', img_prev, shortHelp='previous layer')
        toolbar.AddLabelTool(wx.ID_ABOUT, 'about', img_help, shortHelp='about')
        toolbar.AddLabelTool(wx.ID_EXIT, 'quit', img_quit, shortHelp='quit')
        toolbar.Realize()

        self.Bind(wx.EVT_TOOL, self.OnOpen, id=wx.ID_OPEN)
        self.Bind(wx.EVT_TOOL, self.OnSave, id=wx.ID_SAVE)
        self.Bind(wx.EVT_TOOL, self.OnSlice, id=self.ID_SLICE)
        self.Bind(wx.EVT_TOOL, self.OnNextLayer, id=self.ID_NEXT)
        self.Bind(wx.EVT_TOOL, self.OnPrevLayer, id=self.ID_PREV)
        self.Bind(wx.EVT_TOOL, self.OnAbout, id=wx.ID_ABOUT)
        self.Bind(wx.EVT_TOOL, self.OnQuit, id=wx.ID_EXIT)
        
    def OnNextLayer(self, event):
        if not self.cadmodel.sliced:
            return
        self.cadmodel.next_layer()
        self.left_panel.set_curr_layer(self.cadmodel.curr_layer + 1)
        self.Refresh()

    def OnPrevLayer(self, event):
        if not self.cadmodel.sliced:
            return

        self.cadmodel.prev_layer()
        self.left_panel.set_curr_layer(self.cadmodel.curr_layer + 1)
        self.Refresh()

    def create_panel(self):
        self.left_panel  = ControlPanel(self, self.cadmodel)
        
        self.sp = wx.SplitterWindow(self)
        self.model_panel = wx.Panel(self.sp, style=wx.SUNKEN_BORDER)
        self.path_panel = wx.Panel(self.sp, style=wx.SUNKEN_BORDER)
        self.path_panel.SetBackgroundColour('sky blue')
        self.sp.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGED, self.OnPosChanging)
        
        # Model canvas
        self.model_canvas = ModelCanvas(self.model_panel, self.cadmodel)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.model_canvas, 1, wx.EXPAND)
        self.model_panel.SetSizer(sizer)

        box = wx.BoxSizer(wx.HORIZONTAL)
        box.Add(self.left_panel, 0, wx.EXPAND)
        box.Add(self.sp, 1, wx.EXPAND)
        self.SetSizer(box)

        # Path canvas
        self.path_canvas = PathCanvas(self.path_panel, self.cadmodel)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.path_canvas, 1, wx.EXPAND)
        self.path_panel.SetSizer(sizer)

        self.sp.Initialize(self.model_panel)
        self.sp.SplitVertically(self.model_panel, self.path_panel, 300)
        self.sp.SetMinimumPaneSize(20)

        for i in self.left_panel.dimensionPanel.ids:
            print "YEPPP"
            self.Bind(wx.EVT_TEXT_ENTER, self.OnDimChange, id=self.left_panel.dimensionPanel.ids[i])

    def OnDimChange(self, event):
        self.model_panel.Refresh(False)
        event.Skip()
    
    def OnPosChanging(self, event):
        self.Refresh(False)

    def create_menubar(self):
        menubar = wx.MenuBar()
        for data in self.menu_data():
            label = data[0]
            items = data[1:]
            menubar.Append(self.create_menu(items), label)
        self.SetMenuBar(menubar)    

    def menu_data(self):
        return (("&File", ("&Open\tCtrl+o", "Open CAD file", self.OnOpen, wx.ID_OPEN),
                          ("S&lice\tCtrl+l", "Slice CAD model", self.OnSlice, -1),
                          ("&Save\tCtrl+s", "Save slice result as xml file", self.OnSave, wx.ID_SAVE),  
                          ("", "", "", ""),
                         ("&Quit\tCtrl+q", "Quit", self.OnQuit, wx.ID_EXIT)),
                ("Edit", ("Next Layer\tpgdn", "next layer", self.OnNextLayer, -1),
                         ("Prev Layer\tpgup", "previous layer", self.OnPrevLayer, -1)),
                ("&Help", ("&About", "About this program", self.OnAbout, wx.ID_ABOUT))
                 )
    
    def OnSave(self, event):
        if not self.cadmodel.sliced:
            return

        wildcard = "xml file (*.xml)|*.xml|All files (*.*)|*.*"
        dlg = wx.FileDialog(None, "Save slice data as xml file", os.getcwd(), self.cadname, wildcard, wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            root, ext = os.path.splitext(filename)
            if ext.lower() != '.xml':
                filename = filename + '.xml'
            self.cadmodel.save(filename)
            print 'slicing info is saved in', filename

    def OnAbout(self, event):
        info = wx.AboutDialogInfo()
        info.Name = "Blackcat"
        info.Version = "0.1"
        info.Copyright = "(C) 2009"
        info.Description = "Slice stl CAD model"
        info.Developers = ["Zhigang Liu"]
        info.License = "GPL2"
        wx.AboutBox(info)

    def create_menu(self, menu_data):
        menu = wx.Menu()
        for label, status, handler, id in menu_data:
            if not label:
                menu.AppendSeparator()
                continue
            menu_item = menu.Append(id, label, status)
            self.Bind(wx.EVT_MENU, handler, menu_item)
        return menu

    def OnOpen(self, event):
        wildcard = "CAD std files (*.stl)|*.stl|All files (*.*)|*.*"
        dlg = wx.FileDialog(None, "Open CAD stl file", os.getcwd(), "", wildcard, wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.statusbar.SetStatusText(path)
            print 'open', path
            ok = self.cadmodel.open(path)
            if ok:
                self.model_canvas.create_model()
                self.path_canvas.Refresh()
                self.left_panel.set_dimension(self.cadmodel.dimension)
                basename = os.path.basename(path)
                root, ext = os.path.splitext(basename)
                self.cadname = root
            else:
                wx.MessageBox("Cannot open " + path, 'Error')
        dlg.Destroy()

    def OnSlice(self, event):
        if not self.cadmodel.loaded:
            wx.MessageBox("load a CAD model first", "warning")
            return

        dlg = ParaDialog(self, self.slice_parameter)
        result = dlg.ShowModal()
        if result == wx.ID_OK:
            dlg.get_values()
            print 'slicing...'
            self.cadmodel.queue = Queue.Queue()
            thread.start_new_thread(self.cadmodel.slice, (self.slice_parameter,))
            num_layers = self.cadmodel.queue.get()
            if num_layers > 0:
                pdlg = wx.ProgressDialog("Slicing in progress", "Progress", 
                                         num_layers, 
                                         style=wx.PD_ELAPSED_TIME|
                                               wx.PD_REMAINING_TIME|
                                               wx.PD_AUTO_HIDE|wx.PD_APP_MODAL)
            
                while True:
                    count = self.cadmodel.queue.get()
                    if count == 'done':
                        count = num_layers
                        pdlg.Update(count)
                        break
                    else:
                        pdlg.Update(count)
                pdlg.Destroy()
            
            self.model_canvas.create_model()
            #self.left_panel.set_dimension(self.cadmodel.dimension)
            self.left_panel.set_slice_info(self.slice_parameter)
            self.path_canvas.Refresh()

            if self.cadmodel.sliced:
                self.left_panel.set_num_layer(len(self.cadmodel.layers))
                self.left_panel.set_curr_layer(self.cadmodel.curr_layer + 1)
            else:
                wx.MessageBox("no layers", "Warning")

        dlg.Destroy()

    def OnQuit(self, event):
        self.Close() 

class CharValidator(wx.PyValidator):
    def __init__(self, data, key):
        wx.PyValidator.__init__(self)
        self.Bind(wx.EVT_CHAR, self.OnChar)
        self.data = data
        self.key = key

    def Clone(self):
        return CharValidator(self.data, self.key)
    
    def Validate(self, win):
        text_ctrl = self.GetWindow()
        text = text_ctrl.GetValue()
        if len(text) == 0:
            wx.MessageBox("This field must contain some text!", "Error")
            text_ctrl.SetBackgroundColour('pink')
            text_ctrl.SetFocus()
            text_ctrl.Refresh()
            return False
        else:
            try:
                value = float(text)
            except ValueError:
                wx.MessageBox("must be a number", "Error")  
                text_ctrl.SetBackgroundColour('pink')
                text_ctrl.SetFocus()
                text_ctrl.Refresh()
                return False
            
            #if value <= 0:
            #    wx.MessageBox("value <= 0!", "Error")
            #    text_ctrl.SetBackgroundColour('pink')
            #    text_ctrl.SetFocus()
            #    text_ctrl.Refresh()
            #    return False

        text_ctrl.SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
        text_ctrl.Refresh()
        return True
    
    def TransferToWindow(self):
        text_ctrl = self.GetWindow()
        value = self.data.get(self.key, "")
        text_ctrl.SetValue(value)
        return True

    def TransferFromWindow(self):
        text_ctrl = self.GetWindow()
        self.data[self.key] = text_ctrl.GetValue()
        return True
    
    def OnChar(self, event):
        code = event.GetKeyCode()
        if code < 256: 
            key = chr(code)
            if key in string.letters:
                return
        event.Skip()

class SlicePanel(wx.Panel):
    def __init__(self, parent, data):
        wx.Panel.__init__(self, parent, -1)
        self.data = data
        self.create_controls()

    def create_controls(self):
        #labels = [("Layer height", "1.0", "height"), ("Pitch", "1.0", "pitch"), \
        #          ("Scanning speed", "20", "speed"), ("Fast speed", "20", "fast")]
        labels = [("Layer height (mm)", "1.0", "height"),("Speed (m/s)","0.3","speed"), \
                    ("Starting X (mm)","0.0","global_start_x"),("Starting Y (mm)","0.0","global_start_y"),("Starting Z (mm)","0.0","global_start_z")]
        
        outsizer = wx.BoxSizer(wx.VERTICAL)
        sizer = wx.BoxSizer(wx.VERTICAL)
        outsizer.Add(sizer, 0, wx.ALL, 10)
        box = wx.FlexGridSizer(rows=3, cols=2, hgap=5, vgap=5)
        for label, dvalue, key in labels:
            lbl = wx.StaticText(self, label=label)
            box.Add(lbl, 0, 0)
            txt = wx.TextCtrl(self, -1, dvalue, size=(80, -1), validator=CharValidator(self.data, key))
            box.Add(txt, 0, 0)
        sizer.Add(box, 0, 0)
        
        # slice direction
        lbl = wx.StaticText(self, label="Slice direction")
        box.Add(lbl, 0, 0)

        self.dir_list = ["+X", "-X", "+Y", "-Y", "+Z", "-Z"]
        self.dir_choice = wx.Choice(self, -1, (160, -1), choices=self.dir_list)
        self.dir_choice.SetStringSelection(self.data['direction'])
        box.Add(self.dir_choice, 0, wx.EXPAND)
        
        # scale
        #lbl = wx.StaticText(self, label="Scale factor")
        #box.Add(lbl, 0, 0)
        #scale_txt = wx.TextCtrl(self, -1, "1", size=(80, -1), validator=CharValidator(self.data, "scale"))
        #box.Add(scale_txt, 0, wx.EXPAND)

        self.SetSizer(outsizer)

    def get_direction(self):
        return self.dir_choice.GetStringSelection()

class ParaDialog(wx.Dialog):
    def __init__(self, parent, slice_parameter):
        self.slice_parameter = slice_parameter
        pre = wx.PreDialog()
        pre.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        pre.Create(parent, -1, "Slice parameters")
        self.PostCreate(pre)
        self.create_controls()

    def create_controls(self):
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.panel = SlicePanel(self, self.slice_parameter)
        sizer.Add(self.panel, 0, 0)
        sizer.Add(wx.StaticLine(self), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)
        
        #
        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        btn_sizer.Add((10, 10), 1)
        ok_btn = wx.Button(self, wx.ID_OK)
        ok_btn.SetDefault()
        cancel_btn = wx.Button(self, wx.ID_CANCEL, "Cancel")
        btn_sizer.Add(ok_btn)
        btn_sizer.Add((10,10), 1)
        btn_sizer.Add(cancel_btn)
        btn_sizer.Add((10,10), 1)
        sizer.Add(btn_sizer, 0, wx.EXPAND|wx.ALL, 10)

        self.SetSizer(sizer)
        self.Fit()
    
    def get_values(self):
        self.slice_parameter["direction"] = self.panel.get_direction()

class BlackcatApp(wx.App):
    def __init__(self, redirect=False, filename=None):
        wx.App.__init__(self, redirect, filename)

    def OnInit(self):
        self.frame = BlackcatFrame()
        self.frame.Show()
        return True

if __name__ == '__main__':
    app = BlackcatApp()
    app.MainLoop()