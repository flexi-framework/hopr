#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from model import MainModel
from view import MainView

class App(QtWidgets.QApplication):
    def __init__(self, scriptpath, sys_argv):
       super(App, self).__init__(sys_argv)
       self.model = MainModel()
       self.main_view = MainView(self.model, scriptpath)
       self.main_view.show() #Maximized()
       self.model.gridChanged.emit()

if __name__ == '__main__':
   scriptpath = os.path.dirname(os.path.abspath(sys.argv[0]))

   app = App(scriptpath, sys.argv)
   sys.exit(app.exec_())
