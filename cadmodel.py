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

from error import *
from point import *
from line import *
from functions import *
from facet import *
from layer import *

class CadModel:
    def __init__(self):
        self.init_logger()
        self.loaded = False
        self.curr_layer = -1
        self.sliced = False
        self.dimension = {}
    
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
        self.scale = float(para["scale"])
        self.global_start = Point(float(para["global_start_x"]),float(para["global_start_y"]),float(para["global_start_z"]))
        
        self.scale_model(self.scale)
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
        self.dimension["oldx"] = str(self.xsize)
        self.dimension["oldy"] = str(self.ysize)
        self.dimension["oldz"] = str(self.zsize)
        self.dimension["newx"] = ""
        self.dimension["newy"] = ""
        self.dimension["newz"] = ""

    def set_new_dimension(self):
        self.dimension["newx"] = str(self.xsize)
        self.dimension["newy"] = str(self.ysize)
        self.dimension["newz"] = str(self.zsize)

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