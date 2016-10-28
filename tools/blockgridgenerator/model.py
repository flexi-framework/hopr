# -*- coding: utf-8 -*
import os
import re
import codecs
from PyQt5.QtCore import pyqtSignal, QObject

from block import Block, Blocks
from boundarycondition import BC, BCs

class MainModel(QObject) :
    # Signals to send messages to the GUI
    gridChanged = pyqtSignal()

    def __init__(self) :
        QObject.__init__(self)
        # empty lists (filled with readProdukte and readZettel)
        #self.blocks = []
        self.blocks = Blocks() 
        self.bcs = BCs()
        self.postRefine = 10

    def addBlock(self,xmin,xmax,ymin,ymax,xcells,ycells,bcxmin_name,bcxmax_name,bcymin_name,bcymax_name, row = -1) :
        bcxmin = self.bcs.addBC(bcxmin_name)
        bcxmax = self.bcs.addBC(bcxmax_name)
        bcymin = self.bcs.addBC(bcymin_name)

        bcymax = self.bcs.addBC(bcymax_name)
        if row > -1 :
            self.blocks.item(row).edit(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax)
        else :
            self.blocks.appendRow(Block(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax))
        self.gridChanged.emit()


    def load(self, filename) :
        lines = open(filename, 'r').readlines()
        try :
            self.postRefine = int(lines[0])
            len_bcs, len_blocks = [int(x) for x in lines[1].strip().split()]
            for i in range(len_bcs) :
                l = lines[i+2]
                tmp = l.strip().split()
                name = tmp[0]
                type,state,periodic = [int(x) for x in tmp[1:]]
                self.bcs.addBC(name,type,state,periodic)
            self.blocks.load(lines[2+len_bcs:], self.bcs)
            self.gridChanged.emit()
        except Exception,e :
            print e
            print "Loading of %s failed" % filename

    def save(self, filename, postRefine) :
        self.postRefine = postRefine
        f = open(filename, 'w') 
        f.write("%d\n" % postRefine)
        f.write("%d %d\n" % (self.bcs.rowCount(), self.blocks.rowCount()))
        for bc in self.bcs :
            bc.writeConfig(f)
        for block in self.blocks :
            block.writeConfig(f)
        f.close()


    def export(self, filename, postRefine) :
        self.postRefine = postRefine
        f = open(filename, 'w')

        project = os.path.splitext(os.path.basename(filename))[0]
        f.write("ProjectName  = %s\n" % project)
        f.write("Debugvisu    = T\n")
        f.write("DebugVisuLevel=2\n")
        f.write("NVisu        =1\n")
        f.write("Mode         =1\n")


        f.write("\n")
        nZones = 0
        for block in self.blocks :
            nZones = block.export(f,self.bcs,postRefine, nZones)

        f.write("nZones = %d\n" % nZones)

        for bc in self.bcs :
            f.write("BoundaryName=BC_%s\n" % bc.name.text()) 
            f.write("BoundaryType=(/%d,0,%d,%d/)\n" % (bc.type.value, bc.state.value, bc.periodic.value)) 
        
        f.write("\n")
        f.write("BoundaryName=BC_zminus   \n")
        f.write("BoundaryType=(/1,0,0,1/) \n")
        f.write("BoundaryName=BC_zplus    \n")
        f.write("BoundaryType=(/1,0,0,-1/)\n")
        f.write("vv=(/0.,0.,1./)          \n")

        f.close()



    
    # finish program => close everything
    def close(self) :
        pass

