# -*- coding: utf-8 -*
import os
import re
import codecs
import math
from PyQt5.QtCore import pyqtSignal, QObject

class Block(QObject) :
    def __init__(self, xmin,xmax, ymin,ymax, xcells,ycells, bcxmin=0,bcxmax=0,bcymin=0,bcymax=0) :
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.xcells = xcells
        self.ycells = ycells
        self.bcxmin = bcxmin
        self.bcxmax = bcxmax
        self.bcymin = bcymin
        self.bcymax = bcymax
        self.w = self.xmax - self.xmin
        self.h = self.ymax - self.ymin

    def isValid(self) :
        return self.w > 0 and self.h > 0 and self.xcells > 0 and self.ycells > 0

    def findBottomLeftCornerIndex(self, x, y) :
        i = int(math.floor((x - self.xmin) / (self.w) * self.xcells))
        j = int(math.floor((y - self.ymin) / (self.h) * self.ycells))
        return i, j
    
    def findTopRightCornerIndex(self, x, y) :
        i = int(math.ceil((x - self.xmin) / (self.w) * self.xcells))
        j = int(math.ceil((y - self.ymin) / (self.h) * self.ycells))
        return i, j

    def getCorner(self, i, j) :
        xc = self.xmin + i*(self.xmax-self.xmin)/self.xcells
        yc = self.ymin + j*(self.ymax-self.ymin)/self.ycells
        return xc,yc

    def refine(self, xs,ys, xe,ye) :
        iBL, jBL = self.findBottomLeftCornerIndex(xs, ys)
        xBL, yBL = self.getCorner(iBL,jBL)
        iTR, jTR = self.findTopRightCornerIndex(xe, ye)
        xTR, yTR = self.getCorner(iTR,jTR)

        if not (0 <= iBL <= self.xcells and 0 <= jBL <= self.ycells) :
            print "Anfangspunkt nicht innerhalb des Blocks!" 
            return []
        if not (0 <= iTR <= self.xcells and 0 <= jTR <= self.ycells) :
            print "Endpunkt nicht innerhalb des Blocks!" 
            return []
        

        # generate new block 0-4 
        # self is shrinked to the center block
        # ------------
        # |  3   |   |
        # -------- 2 |
        # |   |  |   |
        # | 0 |------|
        # |   |   1  |
        # ------------
        blocks = []

        # create new blocks with empty BC
        b0 = Block(self.xmin, xBL      , self.ymin, yTR      , iBL            , jTR)
        b1 = Block(xBL      , self.xmax, self.ymin, yBL      , self.xcells-iBL, jBL)
        b2 = Block(xTR      , self.xmax, yBL      , self.ymax, self.xcells-iTR, self.ycells-jBL)
        b3 = Block(self.xmin, xTR      , yTR      , self.ymax, iTR            , self.ycells-jTR)

        # set BCs
        # b0
        b0.bcxmin = self.bcxmin
        b0.bcxmax = 0
        b0.bcymin = self.bcymin
        if (b3.isValid()) : b0.bcymax = 0
        else :              b0.bcymax = self.bcymax
        # b1
        if (b0.isValid()) : b1.bcxmin = 0
        else :              b1.bcxmin = self.bcxmin
        b1.bcxmax = self.bcxmax
        b1.bcymin = self.bcymin
        b1.bcymax = 0
        # b2
        b2.bcxmin = 0
        b2.bcxmax = self.bcxmax
        if (b1.isValid()) : b2.bcymin = 0
        else :              b2.bcymin = self.bcymin
        b2.bcymax = self.bcymax
        # b3
        b3.bcxmin = self.bcxmin
        if (b2.isValid()) : b3.bcxmax = 0
        else :              b3.bcxmax = self.bcxmax
        b3.bcymin = 0
        b3.bcymax = self.bcymax

        if b0.isValid() :
            blocks.append(b0) 
        if b1.isValid() :
            blocks.append(b1) 
        if b2.isValid() :
            blocks.append(b2) 
        if b3.isValid() :
            blocks.append(b3) 

        self.xmin = xBL
        self.xmax = xTR
        self.ymin = yBL
        self.ymax = yTR
        self.xcells = (iTR - iBL)*2
        self.ycells = (jTR - jBL)*2
        self.w = self.xmax - self.xmin
        self.h = self.ymax - self.ymin

        # set new BC of shrinked center block
        if (b0.isValid()) : self.bcxmin = 0
        if (b1.isValid()) : self.bcymin = 0
        if (b2.isValid()) : self.bcxmax = 0
        if (b3.isValid()) : self.bcymax = 0

        return blocks

    def writeConfig(self, f) :
        f.write("%f %f %f %f %d %d %d %d %d %d\n" % (self.xmin, self.xmax, self.ymin, self.ymax, \
                self.xcells, self.ycells, self.bcxmin, self.bcxmax, self.bcymin, self.bcymax))

    def write(self, f, bcs, postRefine) :
        s = "  Corner   =(/"
        s = s + ",".join([str(tmp) for tmp in [self.xmin,self.ymin,0]]) + " ,, "
        s = s + ",".join([str(tmp) for tmp in [self.xmax,self.ymin,0]]) + " ,, "
        s = s + ",".join([str(tmp) for tmp in [self.xmax,self.ymax,0]]) + " ,, "
        s = s + ",".join([str(tmp) for tmp in [self.xmin,self.ymax,0]]) + " ,, "
        s = s + ",".join([str(tmp) for tmp in [self.xmin,self.ymin,1]]) + " ,, "
        s = s + ",".join([str(tmp) for tmp in [self.xmax,self.ymin,1]]) + " ,, "
        s = s + ",".join([str(tmp) for tmp in [self.xmax,self.ymax,1]]) + " ,, "
        s = s + ",".join([str(tmp) for tmp in [self.xmin,self.ymax,1]]) + " /)\n"
        f.write(s)
        f.write("  nElems   =(/" + ",".join([str(tmp) for tmp in [self.xcells*postRefine,self.ycells*postRefine,1]]) + "/)\n")
        f.write("  elemtype =108\n")
        f.write("  BCIndex  =(/" + ",".join([str(tmp) for tmp in [len(bcs)+1,self.bcymin,self.bcxmax,self.bcymax,self.bcxmin,len(bcs)+2]]) + "/)\n")
        f.write("\n")

    def __str__(self) :
        return "[%0.3f,%0.3f] x [%0.3f,%0.3f] : %dx%d" % (self.xmin,self.xmax, self.ymin,self.ymax, self.xcells,self.ycells)

class MainModel(QObject) :
    # Signals to send messages to the GUI
    gridChanged = pyqtSignal()

    def __init__(self) :
        QObject.__init__(self)
        # empty lists (filled with readProdukte and readZettel)
        self.blocks = []
        self.postRefine = 10

    def addBlock(self,xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax) :
        self.blocks.append(Block(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax))
        self.gridChanged.emit()

    def deleteBlock(self, i) :
        del self.blocks[i]
        self.gridChanged.emit()

    def findBlock(self,x,y) :
        for i in range(len(self.blocks)) :
            b = self.blocks[i]
            if b.xmin < x < b.xmax and b.ymin < y < b.ymax :
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


    def getBCs(self) :
        bcs = {}
        for b in self.blocks :
            bcs[b.bcxmin] = 0
            bcs[b.bcxmax] = 0
            bcs[b.bcymin] = 0
            bcs[b.bcymax] = 0

        
        bcs = sorted(bcs.keys())
        try :
            del bcs[bcs.index(0)] # remove inner BCs
        except :
            pass
        return bcs

    def load(self, filename) :
        lines = open(filename, 'r').readlines()
        try :
            self.postRefine = int(lines[0])
            for l in lines[1:] :
                tmp = l.strip().split()
                xmin,xmax,ymin,ymax = [float(s) for s in tmp[0:4]]
                xcells,ycells, bcxmin,bcxmax,bcymin,bcymax = [int(s) for s in tmp[4:]]
                self.blocks.append(Block(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax))
            self.gridChanged.emit()
        except :
            print "Loading of %s failed" % filename

    def save(self, filename, postRefine) :
        self.postRefine = postRefine
        f = open(filename, 'w') 
        f.write("%d\n" % postRefine)
        for b in self.blocks :
            b.writeConfig(f)
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

        bcs = self.getBCs()

        f.write("\n")
        for b in self.blocks :
            b.write(f,bcs,postRefine)

        f.write("nZones = %d\n" % len(self.blocks)) 

        for bc in bcs :
            f.write("BoundaryName=BC_%d\n" % bc) 
            f.write("BoundaryType=(/2,0,0,0/)\n") 
        
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

