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

from cadmodel import *
from error import *
from point import *
from line import *
from functions import *
from facet import *
from layer import *

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

        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnMouseDown)
        self.Bind(wx.EVT_LEFT_UP, self.OnMouseUp)
        self.Bind(wx.EVT_MOTION, self.OnMouseMotion)

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
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.txt_fields = {}
        self.create_controls()

    def create_controls(self):
        box = wx.StaticBox(self, label="Dimension") 
        sizer = wx.StaticBoxSizer(box, wx.HORIZONTAL)
        self.SetSizer(sizer)
        
        label = "Original"
        items = [("X", "oldx"), ("Y", "oldy"), ("Z", "oldz")]
        s1 = self.create_dimension(label, items)
        sizer.Add(s1, 1, wx.EXPAND|wx.ALL, 2)

        label = "Scaled"
        items = [("X", 'newx'), ('Y', 'newy'), ('Z', 'newz')]
        s2 = self.create_dimension(label, items)
        sizer.Add(s2, 1, wx.EXPAND|wx.ALL, 2)

    def create_dimension(self, label, items):
        sizer = wx.BoxSizer(wx.VERTICAL) 
        caption = wx.StaticText(self, label=label)
        sizer.Add(caption, 0, wx.ALIGN_CENTER)

        flex = wx.FlexGridSizer(rows=len(items), cols=2, hgap=2, vgap=2)
        for label, key in items:
            lbl_ctrl = wx.StaticText(self, label=label)
            txt_ctrl = wx.TextCtrl(self, size=(70, -1), style=wx.TE_READONLY)
            flex.Add(lbl_ctrl)
            flex.Add(txt_ctrl, 0, wx.EXPAND)
            self.txt_fields[key] = txt_ctrl
        sizer.Add(flex, 0, wx.EXPAND)
        flex.AddGrowableCol(1, 1)
        return sizer

    def set_values(self, dimension):
        for key in dimension:
            self.txt_fields[key].SetValue(dimension[key])

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
        self.dimensionPanel = DimensionPanel(self)
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
        box2 = wx.BoxSizer(wx.HORIZONTAL)
        axis = ["+X", "-X", "+Y", "-Y", "+Z", "-Z"]
        box2.Add(wx.ComboBox(self, value='+Z', pos=wx.DefaultPosition, size=wx.DefaultSize, choices=axis, style=wx.CB_READONLY), 1)
        sizer.Add(box2, 0, wx.EXPAND)

        self.Bind(wx.EVT_COMBOBOX, self.OnSelect)

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

    def OnSelect(self, event):
        item = event.GetSelection()
        axis = ['+X', '-X', '+Y', '-Y', '+Z', '-Z']
        print "ON SELECT EVENT"
        print axis[item]
        print "THAT's WHAT HAPPENED"
        self.cadmodel.direction = axis[item]
        self.cadmodel.change_direction(axis[item])
        self.cadmodel.calc_dimension()
        self.cadmodel.set_new_dimension()

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
        self.slice_parameter = {"height":"1.0", "speed":"0.3", "pitch": "1.0", "direction":"+Z", "scale":"1", "global_start_x":"0", "global_start_y":"0", "global_start_z":"0"}
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
            self.left_panel.set_dimension(self.cadmodel.dimension)
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
        lbl = wx.StaticText(self, label="Scale factor")
        box.Add(lbl, 0, 0)
        scale_txt = wx.TextCtrl(self, -1, "1", size=(80, -1), validator=CharValidator(self.data, "scale"))
        box.Add(scale_txt, 0, wx.EXPAND)
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