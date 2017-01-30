# -*- coding: utf-8 -*
import math
from PyQt5 import QtCore, QtGui, QtWidgets, uic

class IntItem(QtGui.QStandardItem) :
    def __init__(self, v) :
        QtGui.QStandardItem.__init__(self, str(v))
        self.value = v

    def setData(self, v, role=QtCore.Qt.UserRole + 1) :
        if role == QtCore.Qt.EditRole:
            self.value = int(v)
            return
        return super(IntItem, self).setData(v, role)


    def data(self, role=QtCore.Qt.UserRole + 1):
        if role == QtCore.Qt.DisplayRole:
            return str(self.value) 
        return super(IntItem, self).data(role)

class BC(object) :
    def __init__(self, name, type, state, periodic, color=0) :
        self.name = QtGui.QStandardItem(name)
        self.name.setBackground(QtGui.QBrush(QtCore.Qt.GlobalColor(color + 7)))
        self.type = IntItem(type) 
        self.state = IntItem(state)
        self.periodic = IntItem(periodic)
        self.color = color

    def toRow(self) :
        return [self.name, self.type, self.state, self.periodic]

    def writeConfig(self, f) :
        f.write("%s\n" % str(self))

    def __str__(self) :
        return "%s %d %d %d" %(self.name.text(), self.type.value, self.state.value, self.periodic.value)

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
        self.setHeaderData(0, QtCore.Qt.Horizontal, "Name")
        self.setHeaderData(1, QtCore.Qt.Horizontal, "Type")
        self.setHeaderData(2, QtCore.Qt.Horizontal, "State")
        self.setHeaderData(3, QtCore.Qt.Horizontal, "Periodic")

        self.internalBCs = []

    def addBC(self, name, type=2, state=0, periodic=0) :
        if name == '0' : return None

        for bc in self.internalBCs :
            if bc.name.text() == name :
                return bc

        self.internalBCs.append(BC(name, type, state, periodic, color=len(self.internalBCs)))
        self.appendRow(self.internalBCs[-1].toRow())
        return self.internalBCs[-1]

    def getByName(self, name) :
        for bc in self.internalBCs :
            if bc.name.text() == name :
                return bc
        return None

    def getIndex(self, bc) :
        try :
            return self.internalBCs.index(bc) + 1
        except :
            return 0

    def __iter__(self) :
        return BCsIterator(self)

