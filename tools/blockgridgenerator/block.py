# -*- coding: utf-8 -*
import math
from PyQt5 import QtCore, QtGui, QtWidgets, uic

class Block(QtGui.QStandardItem) :
    def __init__(self, xmin,xmax, ymin,ymax, xcells,ycells, bcxmin=0,bcxmax=0,bcymin=0,bcymax=0) :
        QtGui.QStandardItem.__init__(self)
        self.edit(xmin,xmax, ymin,ymax, xcells,ycells, bcxmin,bcxmax,bcymin,bcymax)

    def edit(self, xmin,xmax, ymin,ymax, xcells,ycells, bcxmin,bcxmax,bcymin,bcymax) :
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

    def isInside(self, x, y) :
        return self.xmin < x < self.xmax and self.ymin < y < self.ymax 

    def drawVerticalLine(self, p, scale, i, penwidth) :
        if i == 0 :
            color = QtCore.Qt.GlobalColor(self.bcxmin + 6)
        elif i == self.xcells :
            color = QtCore.Qt.GlobalColor(self.bcxmax + 6)
        else :
            color = QtCore.Qt.black
        p.setPen(QtGui.QPen(color, penwidth/scale, QtCore.Qt.SolidLine))
        x = self.xmin + self.w/self.xcells * i
        p.drawLine(QtCore.QPointF(x,self.ymin),QtCore.QPointF(x,self.ymax))

    def drawHorizontalLine(self, p, scale, j, penwidth) :
        if j == 0 :
            color = QtCore.Qt.GlobalColor(self.bcymin + 6)
        elif j == self.ycells :
            color = QtCore.Qt.GlobalColor(self.bcymax + 6)
        else :
            color = QtCore.Qt.black
        p.setPen(QtGui.QPen(color, penwidth/scale, QtCore.Qt.SolidLine))
        y = self.ymin + self.h/self.ycells * j
        p.drawLine(QtCore.QPointF(self.xmin,y),QtCore.QPointF(self.xmax,y))

    def draw(self, p, scale, penwidth) :
        print "draw"
        # draw vertical lines
        for i in range(self.xcells+1) :
            self.drawVerticalLine(p, scale, i, penwidth)

        # draw horizontal lines
        for j in range(self.ycells+1) :
            self.drawHorizontalLine(p, scale, j, penwidth)

        for i in range(self.rowCount()) :
            self.child(i).draw(p, scale, penwidth)


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
            self.appendRow(b0) 
        if b1.isValid() :
            self.appendRow(b1) 
        if b2.isValid() :
            self.appendRow(b2) 
        if b3.isValid() :
            self.appendRow(b3) 

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

        #return blocks

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
        try :
            ibcxmin = bcs.index(self.bcxmin) + 1
        except :
            ibcxmin = 0
        try :
            ibcxmax = bcs.index(self.bcxmax) + 1
        except :
            ibcxmax = 0
        try :
            ibcymin = bcs.index(self.bcymin) + 1
        except :
            ibcymin = 0
        try :
            ibcymax = bcs.index(self.bcymax) + 1
        except :
            ibcymax = 0
        f.write("  BCIndex  =(/" + ",".join([str(tmp) for tmp in [len(bcs)+1,ibcymin,ibcxmax,ibcymax,ibcxmin,len(bcs)+2]]) + "/)\n")
        f.write("\n")


    def data(self, role=QtCore.Qt.UserRole + 1):
        if role == QtCore.Qt.DisplayRole:
            return str(self) 
        return super(Block, self).data(role)

    def __str__(self) :
        return "[%0.3f,%0.3f] x [%0.3f,%0.3f] : %dx%d" % (self.xmin,self.xmax, self.ymin,self.ymax, self.xcells,self.ycells)

class BlocksIterator:
    def __init__(self, blocks):
        self.i = 0
        self.blocks = blocks

    def __iter__(self):
        return self

    def next(self):
        if self.i < self.blocks.rowCount():
            i = self.i
            self.i += 1
            return self.blocks.item(i) 
        else:
            raise StopIteration()

class Blocks(QtGui.QStandardItemModel) :
    def __init__(self) :
        QtGui.QStandardItemModel.__init__(self, 0,1)    

    def getTotalSizes(self) :
        xtotalmin = min([block.xmin for block in self])
        ytotalmin = min([block.ymin for block in self])
        xtotalmax = max([block.xmax for block in self])
        ytotalmax = max([block.ymax for block in self])
        return xtotalmin,xtotalmax,ytotalmin,ytotalmax

    def __iter__(self) :
        return BlocksIterator(self)

