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
        self.b_editBlock.clicked.connect(lambda : self.addBlock(edit=True))
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


    def setBCtext(self, lineEdit, bc) :
        if bc :
            lineEdit.setText(bc.name.text())
        else :
            lineEdit.setText("0")

    def selectionChanged_blocks(self) :
        block = self.getSelectedBlock()
        if block :
            self.sb_xmin.setValue(block.xmin)
            self.sb_xmax.setValue(block.xmax)

            self.sb_ymin.setValue(block.ymin)
            self.sb_ymax.setValue(block.ymax)

            self.sb_xcells.setValue(block.xcells)
            self.sb_ycells.setValue(block.ycells)

            self.setBCtext(self.le_bcxmin, block.bcxmin)
            self.setBCtext(self.le_bcxmax, block.bcxmax)

            self.setBCtext(self.le_bcymin, block.bcymin)
            self.setBCtext(self.le_bcymax, block.bcymax)
            self.setEnabledEditButton()

        self.draw()

    def setEnabledEditButton(self) :
        block = self.getSelectedBlock()
        if block :
            self.b_editBlock.setEnabled(not (block.hasChildren() or block.parent()))

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
        if self.model.blocks.rowCount() == 0 : return
        self.startx = e.localPos().x() / self.scale
        h = self.l_draw.height()
        self.starty = (h-e.localPos().y()) / self.scale
        xtotalmin = min([b.xmin for b in self.model.blocks])
        ytotalmin = min([b.ymin for b in self.model.blocks])
        self.startx = self.startx + xtotalmin
        self.starty = self.starty + ytotalmin

    def mouseRelease(self, e) :
        if self.model.blocks.rowCount() == 0 : return
        x = e.localPos().x() / self.scale
        h = self.l_draw.height()
        y = (h-e.localPos().y()) / self.scale
        xtotalmin = min([b.xmin for b in self.model.blocks])
        ytotalmin = min([b.ymin for b in self.model.blocks])
        x = x + xtotalmin
        y = y + ytotalmin
        self.model.blocks.refine(min(x,self.startx), min(y,self.starty), max(x,self.startx), max(y,self.starty))
        self.draw()

    def addBlock(self, edit=False) :
        xmin = self.sb_xmin.value()
        xmax = self.sb_xmax.value()

        ymin = self.sb_ymin.value()
        ymax = self.sb_ymax.value()

        xcells = self.sb_xcells.value()
        ycells = self.sb_ycells.value()

        bcxmin_name = self.le_bcxmin.text().strip().replace(" ","")
        bcxmax_name = self.le_bcxmax.text().strip().replace(" ","")
        bcymin_name = self.le_bcymin.text().strip().replace(" ","")
        bcymax_name = self.le_bcymax.text().strip().replace(" ","")

        if edit :
            index = self.getSelectedIndex(self.tv_blocks)
            if index :
                self.model.addBlock(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin_name,bcxmax_name,bcymin_name,bcymax_name,row=index.row())
        else :
            self.model.addBlock(xmin,xmax,ymin,ymax,xcells,ycells,bcxmin_name,bcxmax_name,bcymin_name,bcymax_name)


    def clearBlocks(self) :
       self.model.blocks.clear()
       self.draw()

    # general routine to get the index of the first selected item of a listview or tableview
    def getSelectedIndex(self, listview) : 
        selection = listview.selectionModel()
        indexes = selection.selectedIndexes()
        if indexes:
            return indexes[0]
        else :
            return None

    def getSelectedBlock(self) :
        index = self.getSelectedIndex(self.tv_blocks)
        if index :
            return self.model.blocks.itemFromIndex(index)
        return None

    def getSelectedBC(self) :
        index = self.getSelectedIndex(self.tv_blocks)
        if index :
            return self.model.blocks.itemFromIndex(index)
        return None

    def keyPressEvent_tv_blocks(self, e) :
        key = e.key()
        if key == QtCore.Qt.Key_Delete :
            # remove selected block
            index = self.getSelectedIndex(self.tv_blocks)
            if index :
                item = self.model.blocks.itemFromIndex(index) 
                parent = item.parent()
                if parent :
                    parent.removeRows(0,parent.rowCount())
                else :
                    self.model.blocks.removeRows(index.row(), 1)
                self.draw()

    def draw(self) :
        self.setEnabledEditButton()
        self.model.bcs.sort(0,0)
        self.sb_postRefine.setValue(self.model.postRefine)
      
        w = self.l_draw.width()-2
        h = self.l_draw.height()-2

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
            block.draw(p, self.scale, selected_block)
        p.end()

        self.l_draw.setPicture(pic) 
        self.l_draw.show()

    # end of programm (X clicked)
    def closeEvent(self, event) :
        self.model.close()
