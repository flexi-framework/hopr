# -*- coding: utf-8 -*
import math
from PyQt5 import QtCore, QtGui, QtWidgets, uic

class IntItem(QtGui.QStandardItem) :
    def __init__(self, v) :
        QtGui.QStandardItem.__init__(self, str(v))
        self.value = v

    def setData(self, v, role=QtCore.Qt.UserRole + 1) :
        self.value = int(v)

    def data(self, role=QtCore.Qt.UserRole + 1):
        if role == QtCore.Qt.DisplayRole:
            return str(self.value) 
        return super(IntItem, self).data(role)

class BC(object) :
    def __init__(self, no, type=2, state=0, periodic=0) :
        self.no = IntItem(no)
        self.type = IntItem(type) 
        self.state = IntItem(state)
        self.periodic = IntItem(periodic)

    def toRow(self) :
        return [self.no, self.type, self.state, self.periodic]

    def writeConfig(self, f) :
        f.write("%d %d %d %d\n" %(self.no, self.type, self.state, self.periodic))

    def __str__(self) :
        return "%d %d %d %d" %(self.no.value, self.type.value, self.state.value, self.periodic.value)

class BCsIterator:
    def __init__(self, bcs):
        self.i = 0
        self.bcs = bcs

    def __iter__(self):
        return self

    def next(self):
        if self.i < self.bcs.rowCount():
            i = self.i
            self.i += 1
            return self.bcs.internalBCs[i] 
        else:
            raise StopIteration()

class BCs(QtGui.QStandardItemModel) :
    def __init__(self) :
        QtGui.QStandardItemModel.__init__(self, 0,4)    
        # set the headers of the columns
        self.setHeaderData(0, QtCore.Qt.Horizontal, "Nummer")
        self.setHeaderData(1, QtCore.Qt.Horizontal, "Type")
        self.setHeaderData(2, QtCore.Qt.Horizontal, "State")
        self.setHeaderData(3, QtCore.Qt.Horizontal, "Periodic")

        self.internalBCs = []

    def addBC(self, bc_no) :
        for bc in self.internalBCs :
            if bc.no == bc_no :
                return

        self.internalBCs.append(BC(bc_no))
        self.appendRow(self.internalBCs[-1].toRow())


    def __iter__(self) :
        return BCsIterator(self)

