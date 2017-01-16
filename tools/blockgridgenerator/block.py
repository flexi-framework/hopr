# -*- coding: utf-8 -*
import math
from PyQt5 import QtCore, QtGui, QtWidgets, uic

class Block(QtGui.QStandardItem) :
    def __init__(self, xmin,xmax, ymin,ymax, xcells,ycells, bcxmin=None,bcxmax=None,bcymin=None,bcymax=None) :
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

    def getColor(self,i=None,j=None) :
        if 0 < i < self.xcells : return QtCore.Qt.black
        if 0 < j < self.ycells : return QtCore.Qt.black
        if i == 0 and self.bcxmin :
            return QtCore.Qt.GlobalColor(self.bcxmin.color + 7)
        if i == self.xcells and self.bcxmax :
            return QtCore.Qt.GlobalColor(self.bcxmax.color + 7)
        if j == 0 and self.bcymin :
            return QtCore.Qt.GlobalColor(self.bcymin.color + 7)
        if j == self.ycells and self.bcymax :
            return QtCore.Qt.GlobalColor(self.bcymax.color + 7)
        return  QtCore.Qt.lightGray

    def drawVerticalLine(self, p, scale, i, penwidth) :
        color = self.getColor(i=i)
        p.setPen(QtGui.QPen(color, penwidth/scale, QtCore.Qt.SolidLine))
        x = self.xmin + self.w/self.xcells * i
        p.drawLine(QtCore.QPointF(x,self.ymin),QtCore.QPointF(x,self.ymax))

    def drawHorizontalLine(self, p, scale, j, penwidth) :
        color = self.getColor(j=j)
        p.setPen(QtGui.QPen(color, penwidth/scale, QtCore.Qt.SolidLine))
        y = self.ymin + self.h/self.ycells * j
        p.drawLine(QtCore.QPointF(self.xmin,y),QtCore.QPointF(self.xmax,y))

    def draw(self, p, scale, selected_block, penwidth = 1.) :
        if self == selected_block :
            penwidth = 2.
        if self.hasChildren() :
            for i in range(self.rowCount()) :
                self.child(i).draw(p, scale, selected_block, penwidth)
        else :
            # draw vertical lines
            for i in range(self.xcells+1) :
                self.drawVerticalLine(p, scale, i, penwidth)

            # draw horizontal lines
            for j in range(self.ycells+1) :
                self.drawHorizontalLine(p, scale, j, penwidth)


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
        # ------------
        # |  3   |   |
        # -------- 2 |
        # |   |  |   |
        # | 0 |------|
        # |   |   1  |
        # ------------

        # create new blocks with inner BCs
        b0 = Block(self.xmin, xBL      , self.ymin, yTR      , iBL            , jTR)
        b1 = Block(xBL      , self.xmax, self.ymin, yBL      , self.xcells-iBL, jBL)
        b2 = Block(xTR      , self.xmax, yBL      , self.ymax, self.xcells-iTR, self.ycells-jBL)
        b3 = Block(self.xmin, xTR      , yTR      , self.ymax, iTR            , self.ycells-jTR)

        # set BCs
        # b0
        b0.bcxmin = self.bcxmin
        b0.bcxmax = None
        b0.bcymin = self.bcymin
        if (b3.isValid()) : b0.bcymax = None
        else :              b0.bcymax = self.bcymax

        # b1
        if (b0.isValid()) : b1.bcxmin = None
        else :              b1.bcxmin = self.bcxmin
        b1.bcxmax = self.bcxmax
        b1.bcymin = self.bcymin
        b1.bcymax = None

        # b2
        b2.bcxmin = None
        b2.bcxmax = self.bcxmax
        if (b1.isValid()) : b2.bcymin = None
        else :              b2.bcymin = self.bcymin
        b2.bcymax = self.bcymax

        # b3
        b3.bcxmin = self.bcxmin
        if (b2.isValid()) : b3.bcxmax = None
        else :              b3.bcxmax = self.bcxmax
        b3.bcymin = None
        b3.bcymax = self.bcymax

        # Append valid blocks (non-empty)
        if b0.isValid() :
            self.appendRow(b0) 
        if b1.isValid() :
            self.appendRow(b1) 
        if b2.isValid() :
            self.appendRow(b2) 
        if b3.isValid() :
            self.appendRow(b3) 

        # generate center block with inner BC
        bcenter = Block(xBL,xTR,yBL,yTR, (iTR - iBL)*2,(jTR - jBL)*2)

        # set BC of center block, if the adjacent block is not valid
        if not b0.isValid() : bcenter.bcxmin = self.bcxmin
        if not b1.isValid() : bcenter.bcymin = self.bcymin
        if not b2.isValid() : bcenter.bcxmax = self.bcxmax
        if not b3.isValid() : bcenter.bcymax = self.bcymax

        # Append center block
        self.appendRow(bcenter)

    def getBCtext(self, bc) :
        if bc :
            return bc.name.text()
        else :
            return "0"

    def writeConfig(self, f, indent = "") :
        bcxmin_name = self.getBCtext(self.bcxmin)
        bcxmax_name = self.getBCtext(self.bcxmax)
        bcymin_name = self.getBCtext(self.bcymin)
        bcymax_name = self.getBCtext(self.bcymax)
        f.write("%s%f %f %f %f %d %d %s %s %s %s\n" % (indent, self.xmin, self.xmax, self.ymin, self.ymax, \
                self.xcells, self.ycells, bcxmin_name, bcxmax_name, bcymin_name, bcymax_name))
        if self.hasChildren() :
            for i in range(self.rowCount()) :
                self.child(i).writeConfig(f, indent + "   ")

    def export(self, f, bcs, postRefine, nZones) :
        if self.hasChildren() :
            for i in range(self.rowCount()) :
                nZones = self.child(i).export(f, bcs, postRefine, nZones)
            return nZones
        else :
            s = "Corner   =(/"
            s = s + ",".join([str(tmp) for tmp in [self.xmin,self.ymin,0]]) + " ,, "
            s = s + ",".join([str(tmp) for tmp in [self.xmax,self.ymin,0]]) + " ,, "
            s = s + ",".join([str(tmp) for tmp in [self.xmax,self.ymax,0]]) + " ,, "
            s = s + ",".join([str(tmp) for tmp in [self.xmin,self.ymax,0]]) + " ,, "
            s = s + ",".join([str(tmp) for tmp in [self.xmin,self.ymin,1]]) + " ,, "
            s = s + ",".join([str(tmp) for tmp in [self.xmax,self.ymin,1]]) + " ,, "
            s = s + ",".join([str(tmp) for tmp in [self.xmax,self.ymax,1]]) + " ,, "
            s = s + ",".join([str(tmp) for tmp in [self.xmin,self.ymax,1]]) + " /)\n"
            f.write(s)
            f.write("nElems   =(/" + ",".join([str(tmp) for tmp in [self.xcells*postRefine,self.ycells*postRefine,1]]) + "/)\n")
            f.write("elemtype =108\n")
            ibcxmin = bcs.getIndex(self.bcxmin)
            ibcxmax = bcs.getIndex(self.bcxmax)
            ibcymin = bcs.getIndex(self.bcymin)
            ibcymax = bcs.getIndex(self.bcymax)
            ibczmin = len(bcs.internalBCs) + 1
            ibczmax = len(bcs.internalBCs) + 2
            f.write("  BCIndex  =(/" + ",".join([str(tmp) for tmp in [ibczmin,ibcymin,ibcxmax,ibcymax,ibcxmin,ibczmax]]) + "/)\n")
            f.write("\n")
            return nZones + 1


    def getBlockOfXY(self, x,y) :
        if self.hasChildren() :
            for i in range(self.rowCount()) :
                block = self.child(i).getBlockOfXY(x,y)
                if block : return block
            return None
        if self.isInside(x,y) :
            return self
        else :
            return None

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
   
def parseLine(line) :
    for indent in range(len(line)) :
        if line[indent] != " ": break
    tmp = line.strip().split()
    coords = [float(x) for x in tmp[0:4]]
    cells = [int(x) for x in tmp[4:6]]
    bcnames = tmp[6:]
    return indent/3, coords, cells, bcnames

class Blocks(QtGui.QStandardItemModel) :
    def __init__(self) :
        QtGui.QStandardItemModel.__init__(self, 0,1)    
        self.setHeaderData(0, QtCore.Qt.Horizontal, "Blocks")

    def getTotalSizes(self) :
        xtotalmin = min([block.xmin for block in self])
        ytotalmin = min([block.ymin for block in self])
        xtotalmax = max([block.xmax for block in self])
        ytotalmax = max([block.ymax for block in self])
        return xtotalmin,xtotalmax,ytotalmin,ytotalmax

    def refine(self, xs,ys, xe,ye) :
        self.setHeaderData(0, QtCore.Qt.Horizontal, "Blocks")
        block_s = None
        block_e = None
        # get block of start point
        for block in self :
            tmp = block.getBlockOfXY(xs,ys)
            if tmp : 
                block_s = tmp
                break
        # get block of end point
        for block in self :
            tmp = block.getBlockOfXY(xe,ye)
            if tmp : 
                block_e = tmp
                break

        if not block_s == block_e :
            print "Anfangs- und Endpunkt nicht im gleichen Block"
        else :
            if block_s :
                block_s.refine(xs,ys,xe,ye)

    def loadBlock(self, bcs, block, coords, cells, bcnames) :
        bcxmin = bcs.getByName(bcnames[0])
        bcxmax = bcs.getByName(bcnames[1])
        bcymin = bcs.getByName(bcnames[2])
        bcymax = bcs.getByName(bcnames[3])

        newBlock = Block(coords[0], coords[1], coords[2], coords[3], cells[0], cells[1], bcxmin, bcxmax, bcymin, bcymax)
        if block : 
            block.appendRow(newBlock)
        else :
            self.appendRow(newBlock)
        return newBlock

    def load(self, lines, bcs) :
        i = 0
        lastindent = 0
        block = None
        for line in lines :
            indent,coords,cells,bcnames = parseLine(line)
            if indent <= lastindent :
                if block :
                    block = block.parent()
            if indent < lastindent :
                if block :
                    block = block.parent()
            block = self.loadBlock(bcs, block, coords, cells, bcnames)
            lastindent = indent

    def __iter__(self) :
        return BlocksIterator(self)

