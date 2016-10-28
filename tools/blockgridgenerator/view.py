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

        #self.tm_blocks = QtGui.QStandardItemModel(0,1)
        self.tv_blocks.setModel(self.model.blocks)
        self.tv_blocks.keyPressEvent = self.keyPressEvent_tv_blocks
        self.tv_blocks.selectionModel().selectionChanged.connect(self.selectionChanged_blocks)

        # do not show vertical headers
        self.tv_bcs.setModel(self.model.bcs)
        self.tv_bcs.verticalHeader().setVisible(False)
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

    def selectionChanged_blocks(self) :
        block = self.getSelectedBlock()
        if block :
            self.sb_xmin.setValue(block.xmin)
            self.sb_xmax.setValue(block.xmax)

            self.sb_ymin.setValue(block.ymin)
            self.sb_ymax.setValue(block.ymax)

            self.sb_xcells.setValue(block.xcells)
            self.sb_ycells.setValue(block.ycells)

            self.sb_bcxmin.setValue(block.bcxmin)
            self.sb_bcxmax.setValue(block.bcxmax)

            self.sb_bcymin.setValue(block.bcymin)
            self.sb_bcymax.setValue(block.bcymax)

        self.draw()

    def selectionChanged_bcs(self) :
        self.draw()

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
        index = self.getSelectedItem(self.tv_blocks)
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

            self.model.addBlock(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin,bcxmax,bcymin,bcymax,edit=index.row())

    def clearBlocks(self) :
       self.model.blocks.clear()
       self.draw()

    # general routine to get the index of the first selected item of a listview or tableview
    def getSelectedItem(self, listview) : 
        selection = listview.selectionModel()
        indexes = selection.selectedIndexes()
        if indexes:
            return indexes[0]
        else :
            return None

    def getSelectedBlock(self) :
        index = self.getSelectedItem(self.tv_blocks)
        if index :
            return self.model.blocks.itemFromIndex(index)
        return None

    def getSelectedBC(self) :
        index = self.getSelectedItem(self.tv_blocks)
        if index :
            return self.model.blocks.itemFromIndex(index)
        return None

    def keyPressEvent_tv_blocks(self, e) :
        key = e.key()
        if key == QtCore.Qt.Key_Delete :
            # remove selected block
            index = self.getSelectedItem(self.tv_blocks)
            if index :
                self.model.deleteBlock(index.row())

    def draw(self) :
        self.sb_postRefine.setValue(self.model.postRefine)
      
        w = self.l_draw.width()-1
        h = self.l_draw.height()-1

        if self.model.blocks.rowCount() > 0 :
            totalsizes = self.model.blocks.getTotalSizes()
            wpic = totalsizes[1] - totalsizes[0]
            hpic = totalsizes[3] - totalsizes[2]

            self.scale = min(w/wpic, h/hpic)
        else :
            self.scale = 1

        pic = QtGui.QPicture()
        p = QtGui.QPainter()
        p.begin(pic)
        p.scale(self.scale,-self.scale)

        selected_block = self.getSelectedBlock()
        selected_bc = self.getSelectedBC()

        for block in self.model.blocks :
            penwidth = 1.
            if block == selected_block :
                penwidth = 2.
            block.draw(p, self.scale, penwidth)
        p.end()

        self.l_draw.setPicture(pic) 
        self.l_draw.show()

    # end of programm (X clicked)
    def closeEvent(self, event) :
        self.model.close()
