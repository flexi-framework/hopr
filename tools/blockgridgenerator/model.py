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

    def addBlock(self,xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax, row = -1) :
        if row > -1 :
            self.blocks.item(i).edit(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax)
        else :
            self.blocks.appendRow(Block(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax))
        self.bcs.addBC(bcxmin)
        self.bcs.addBC(bcxmax)
        self.bcs.addBC(bcymin)
        self.bcs.addBC(bcymax)
        self.gridChanged.emit()

    def deleteBlock(self, i) :
        self.blocks.removeRows(i, 1)
        self.gridChanged.emit()

    def findBlock(self,x,y) :
        for block in self.blocks : 
            if block.isInside(x,y) :
                return block
        return -1

    def refineBlock(self,xs,ys,xe,ye) :
        block_s = self.findBlock(xs,ys)
        block_e = self.findBlock(xe,ye)
        if block_s != block_e :
            print "Anfangs- und Endpunkt nicht im gleichen Block"
            return

        block_s.refine(xs,ys,xe,ye)
        self.gridChanged.emit()

    def load(self, filename) :
        lines = open(filename, 'r').readlines()
        try :
            self.postRefine = int(lines[0])
            len_blocks, len_bcs = [int(x) for x in lines[1].strip().split()]
            print len_blocks, len_bcs
            for i in range(len_blocks) :
                l = lines[i+2]
                tmp = l.strip().split()
                xmin,xmax,ymin,ymax = [float(s) for s in tmp[0:4]]
                xcells,ycells, bcxmin,bcxmax,bcymin,bcymax = [int(s) for s in tmp[4:]]
                self.blocks.appendRow(Block(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax))
            for i in range(len_bcs) :
                l = lines[i+2+len_blocks]
                no,type,state,periodic = [int(x) for x in l.strip().split()]
                print no, type, state, periodic
                self.bcs[no] = BC(no,type,state,periodic)
            self.gridChanged.emit()
        except Exception,e :
            print e
            print "Loading of %s failed" % filename

    def save(self, filename, postRefine) :
        self.postRefine = postRefine
        f = open(filename, 'w') 
        f.write("%d\n" % postRefine)
        f.write("%d %d\n" % (self.blocks.rowCount(), len(self.bcs)))
        for i in range(self.blocks.rowCount()) :
            self.blocks.item(i).writeConfig(f)
        for bc_no in sorted(self.bcs.keys()) :
            self.bcs[bc_no].writeConfig(f)
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
        for i in range(self.blocks.rowCount()) :
            self.blocks.item(i).write(f,sorted(self.bcs.keys()),postRefine)

        f.write("nZones = %d\n" % len(self.blocks)) 

        for bc_no in sorted(self.bcs.keys()) :
            bc = self.bcs[bc_no]
            f.write("BoundaryName=BC_%d\n" % bc.no) 
            f.write("BoundaryType=(/%d,0,%d,%d/)\n" % (bc.type, bc.state, bc.periodic)) 
        
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

