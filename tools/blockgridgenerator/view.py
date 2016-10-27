# -*- coding: utf-8 -*
import sys
import os
from PyQt5 import QtCore, QtGui, QtWidgets, uic

# Main View
class MainView(QtWidgets.QMainWindow):
    def __init__(self, model, scriptpath):
        QtWidgets.QMainWindow.__init__(self, None)
        uic.loadUi(os.path.join(scriptpath,'mainwindow.ui'), self)
        self.model = model
        # connect the signals from the model with functions
        self.model.gridChanged.connect(self.draw)
        self.scale = 1

        self.b_addBlock.clicked.connect(self.addBlock)
        self.b_editBlock.clicked.connect(self.editBlock)
        self.b_clearBlocks.clicked.connect(self.clearBlocks)
        self.b_load.clicked.connect(self.load)
        self.b_save.clicked.connect(self.save)
        self.b_export.clicked.connect(self.export)

        self.tm_blocks = QtGui.QStandardItemModel(0,1)
        self.lv_blocks.setModel(self.tm_blocks)
        self.lv_blocks.keyPressEvent = self.keyPressEvent_lv_blocks
        self.lv_blocks.selectionModel().selectionChanged.connect(self.selectionChanged_blocks)

        self.tm_bcs = QtGui.QStandardItemModel(0,4)
        # set the headers of the columns
        self.tm_bcs.setHeaderData(0, QtCore.Qt.Horizontal, "Nummer")
        self.tm_bcs.setHeaderData(1, QtCore.Qt.Horizontal, "Type")
        self.tm_bcs.setHeaderData(2, QtCore.Qt.Horizontal, "State")
        self.tm_bcs.setHeaderData(3, QtCore.Qt.Horizontal, "Periodic")
        self.tm_bcs.dataChanged.connect(self.bcsChanged)
        # do not show vertical headers
        self.tv_bcs.verticalHeader().setVisible(False)
        self.tv_bcs.setModel(self.tm_bcs)
        self.tv_bcs.setColumnWidth(0,70)
        self.tv_bcs.setColumnWidth(1,70)
        self.tv_bcs.setColumnWidth(2,70)
        self.tv_bcs.setColumnWidth(3,70)
        self.tv_bcs.selectionModel().selectionChanged.connect(self.selectionChanged_bcs)



        self.l_draw.mousePressEvent = self.mousePress
        self.l_draw.mouseReleaseEvent = self.mouseRelease
        self.l_draw.setMargin(1)

        self.loadname = ""
        self.exportname = ""

    def bcsChanged(self, topLeft, bottomRight, parent) :
        row    = topLeft.row()
        column = topLeft.column()
        
        try :
            selected_bc = int(self.tm_bcs.item(row,0).text())
            # get the new value
            newValue = int(self.tm_bcs.data(topLeft))
        except :
            self.draw()
            return

        if column == 0 :
            # undo change of no
            self.redraw(selected_bc = selected_bc)
        elif column == 1 :
            self.model.changeBCType(selected_bc, newValue) 
        elif column == 2 :
            self.model.changeBCState(selected_bc, newValue) 
        elif column == 3 :
            self.model.changeBCPeriodic(selected_bc, newValue) 

    def selectionChanged_blocks(self) :
        index = self.getSelectedItem(self.lv_blocks)
        if index :
            b = self.model.blocks[index.row()]
            
            self.sb_xmin.setValue(b.xmin)
            self.sb_xmax.setValue(b.xmax)

            self.sb_ymin.setValue(b.ymin)
            self.sb_ymax.setValue(b.ymax)

            self.sb_xcells.setValue(b.xcells)
            self.sb_ycells.setValue(b.ycells)

            self.sb_bcxmin.setValue(b.bcxmin)
            self.sb_bcxmax.setValue(b.bcxmax)

            self.sb_bcymin.setValue(b.bcymin)
            self.sb_bcymax.setValue(b.bcymax)

            self.draw(selected_block = b)

    def selectionChanged_bcs(self) :
        index = self.getSelectedItem(self.tv_bcs)
        if index :
            selected_bc = int(self.tm_bcs.item(index.row(),0).text())
            self.draw(selected_bc = selected_bc)

    def load(self) :
        filename = QtWidgets.QFileDialog.getOpenFileName(self, "Load file", self.loadname, "*.config");
        if os.path.splitext(filename[0])[1] :
            self.loadname = filename[0] 
        else :
            self.loadname = os.path.splitext(filename[0])[0] + os.path.splitext(filename[1])[1]

        if self.loadname :
            self.model.load(self.loadname)

    def save(self) :
        filename = QtWidgets.QFileDialog.getSaveFileName(self, "Save file", self.loadname, "*.config");
        if os.path.splitext(filename[0])[1] :
            self.loadname = filename[0] 
        else :
            self.loadname = os.path.splitext(filename[0])[0] + os.path.splitext(filename[1])[1]

        postRefine = self.sb_postRefine.value()
        if self.loadname :
            self.model.save(self.loadname,postRefine)

    def export(self) :
        filename = QtWidgets.QFileDialog.getSaveFileName(self, "Export file", self.exportname, "*.ini");
        if os.path.splitext(filename[0])[1] :
            self.exportname = filename[0] 
        else :
            self.exportname = os.path.splitext(filename[0])[0] + os.path.splitext(filename[1])[1]
        postRefine = self.sb_postRefine.value()

        if self.exportname :
            self.model.export(self.exportname,postRefine)


    def mousePress(self, e) :
        self.startx = e.localPos().x() / self.scale
        h = self.l_draw.height()
        self.starty = (h-e.localPos().y()) / self.scale
        xtotalmin = min([b.xmin for b in self.model.blocks])
        ytotalmin = min([b.ymin for b in self.model.blocks])
        self.startx = self.startx + xtotalmin
        self.starty = self.starty + ytotalmin

    def mouseRelease(self, e) :
        x = e.localPos().x() / self.scale
        h = self.l_draw.height()
        y = (h-e.localPos().y()) / self.scale
        xtotalmin = min([b.xmin for b in self.model.blocks])
        ytotalmin = min([b.ymin for b in self.model.blocks])
        x = x + xtotalmin
        y = y + ytotalmin
        self.model.refineBlock(min(x,self.startx), min(y,self.starty), max(x,self.startx), max(y,self.starty))

    def addBlock(self) :
        xmin = self.sb_xmin.value()
        xmax = self.sb_xmax.value()

        ymin = self.sb_ymin.value()
        ymax = self.sb_ymax.value()

        xcells = self.sb_xcells.value()
        ycells = self.sb_ycells.value()

        bcxmin = self.sb_bcxmin.value()
        bcxmax = self.sb_bcxmax.value()

        bcymin = self.sb_bcymin.value()
        bcymax = self.sb_bcymax.value()

        self.model.addBlock(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax)

    def editBlock(self) :
        index = self.getSelectedItem(self.lv_blocks)
        if index :

            xmin = self.sb_xmin.value()
            xmax = self.sb_xmax.value()

            ymin = self.sb_ymin.value()
            ymax = self.sb_ymax.value()

            xcells = self.sb_xcells.value()
            ycells = self.sb_ycells.value()

            bcxmin = self.sb_bcxmin.value()
            bcxmax = self.sb_bcxmax.value()

            bcymin = self.sb_bcymin.value()
            bcymax = self.sb_bcymax.value()

            self.model.editBlock(index.row(), xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax)

    def clearBlocks(self) :
       self.model.blocks = []
       self.draw()

    # general routine to get the index of the first selected item of a listview or tableview
    def getSelectedItem(self, listview) : 
        selection = listview.selectionModel()
        indexes = selection.selectedIndexes()
        if indexes:
            return indexes[0]
        else :
            return None

    def keyPressEvent_lv_blocks(self, e) :
        key = e.key()
        if key == QtCore.Qt.Key_Delete :
            # remove selected block
            index = self.getSelectedItem(self.lv_blocks)
            if index :
                self.model.deleteBlock(index.row())


    def draw(self, selected_block=None, selected_bc=None) :
        if not selected_block :
            self.sb_postRefine.setValue(self.model.postRefine)
            self.tm_blocks.removeRows(0, self.tm_blocks.rowCount())
            for block in self.model.blocks :
                self.tm_blocks.appendRow([QtGui.QStandardItem(str(block))])

        if not selected_bc :
            self.tm_bcs.removeRows(0, self.tm_bcs.rowCount())
            for no in sorted(self.model.bcs.keys()) :
                bc = self.model.bcs[no]
                row = [str(bc.no), str(bc.type),str(bc.state), str(bc.periodic)]
                row = [QtGui.QStandardItem(x) for x in row]
                self.tm_bcs.appendRow(row)

        w = self.l_draw.width()-1
        h = self.l_draw.height()-1

        if len(self.model.blocks) > 0 :
            xtotalmin = min([block.xmin for block in self.model.blocks])
            ytotalmin = min([block.ymin for block in self.model.blocks])
            xtotalmax = max([block.xmax for block in self.model.blocks])
            ytotalmax = max([block.ymax for block in self.model.blocks])
            wpic = xtotalmax - xtotalmin
            hpic = ytotalmax - ytotalmin

            self.scale = min(w/wpic, h/hpic)
        else :
            self.scale = 1

        pic = QtGui.QPicture()
        p = QtGui.QPainter()
        p.begin(pic)
        p.scale(self.scale,-self.scale)

        for block in self.model.blocks :
            penwidth = 1.
            if block == selected_block :
                penwidth = 2.
            for i in range(block.xcells+1) :
                x = block.xmin + block.w/block.xcells * i
                if i == 0 :
                    bc = block.bcxmin
                    bold = 1.
                    if bc == selected_bc : bold = 3.
                    p.setPen(QtGui.QPen(QtCore.Qt.GlobalColor(bc + 6), bold*penwidth/self.scale, QtCore.Qt.SolidLine))
                elif i == block.xcells :
                    bc = block.bcxmax
                    bold = 1.
                    if bc == selected_bc : bold = 3.
                    p.setPen(QtGui.QPen(QtCore.Qt.GlobalColor(bc + 6), bold*penwidth/self.scale, QtCore.Qt.SolidLine))
                else :
                    p.setPen(QtGui.QPen(QtCore.Qt.black, penwidth/self.scale, QtCore.Qt.SolidLine))

                p.drawLine(QtCore.QPointF(x,block.ymin),QtCore.QPointF(x,block.ymax))
            for j in range(block.ycells+1) :
                if j == 0 :
                    bc = block.bcymin
                    bold = 1.
                    if bc == selected_bc : bold = 3.
                    p.setPen(QtGui.QPen(QtCore.Qt.GlobalColor(bc + 6), bold*penwidth/self.scale, QtCore.Qt.SolidLine))
                elif j == block.ycells :
                    bc = block.bcymax
                    bold = 1.
                    if bc == selected_bc : bold = 3.
                    p.setPen(QtGui.QPen(QtCore.Qt.GlobalColor(bc + 6), bold*penwidth/self.scale, QtCore.Qt.SolidLine))
                else :
                    p.setPen(QtGui.QPen(QtCore.Qt.black, penwidth/self.scale, QtCore.Qt.SolidLine))

                y = block.ymin + block.h/block.ycells * j
                p.drawLine(QtCore.QPointF(block.xmin,y),QtCore.QPointF(block.xmax,y))

        p.end()

        self.l_draw.setPicture(pic) 
        self.l_draw.show()

    # end of programm (X clicked)
    def closeEvent(self, event) :
        self.model.close()
