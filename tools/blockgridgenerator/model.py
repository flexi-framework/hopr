# -*- coding: utf-8 -*
import os
import re
import codecs
from PyQt5.QtCore import pyqtSignal, QObject

from block import Block

class BC() :
    def __init__(self, no, type=2, state=0, periodic=0) :
        self.no = no
        self.type = type
        self.state = state
        self.periodic = periodic

    def writeConfig(self, f) :
        print "write bc", self.no
        f.write("%d %d %d %d\n" %(self.no, self.type, self.state, self.periodic))

class MainModel(QObject) :
    # Signals to send messages to the GUI
    gridChanged = pyqtSignal()

    def __init__(self) :
        QObject.__init__(self)
        # empty lists (filled with readProdukte and readZettel)
        self.blocks = []
        self.bcs = {}
        self.postRefine = 10

    def addBlock(self,xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax) :
        self.blocks.append(Block(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax))
        self.addBC(bcxmin)
        self.addBC(bcxmax)
        self.addBC(bcymin)
        self.addBC(bcymax)
        self.gridChanged.emit()

    def editBlock(self,i,xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax) :
        del self.blocks[i]
        block = Block(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax)
        self.addBC(bcxmin)
        self.addBC(bcxmax)
        self.addBC(bcymin)
        self.addBC(bcymax)
        self.blocks.insert(i, block)
        self.gridChanged.emit()

    def deleteBlock(self, i) :
        del self.blocks[i]
        self.gridChanged.emit()

    def addBC(self, bc_no) :
        if self.bcs.has_key(bc_no) : return # bc already present
        else :
            self.bcs[bc_no] = BC(bc_no)

    def findBlock(self,x,y) :
        for i in range(len(self.blocks)) :
            block = self.blocks[i]
            if block.xmin < x < block.xmax and block.ymin < y < block.ymax :
                return i
        return -1

    def refineBlock(self,xs,ys,xe,ye) :
        i  = self.findBlock(xs,ys)
        ie = self.findBlock(xe,ye)
        if i != ie :
            print "Anfangs- und Endpunkt nicht im gleichen Block"
            return

        newblocks = self.blocks[i].refine(xs,ys,xe,ye)
        a = 0
        for nb in newblocks :
            a = a+1
            self.blocks.insert(i+a, nb)

        self.gridChanged.emit()

    def changeBCType(self, bc_no, v) :
        self.bcs[bc_no].type = v

    def changeBCState(self, bc_no, v) :
        self.bcs[bc_no].state = v

    def changeBCPeriodic(self, bc_no, v) :
        self.bcs[bc_no].periodic = v

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
                self.blocks.append(Block(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax))
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
        f.write("%d %d\n" % (len(self.blocks), len(self.bcs)))
        for block in self.blocks :
            block.writeConfig(f)
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
        for block in self.blocks :
            block.write(f,sorted(self.bcs.keys()),postRefine)

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

