import sys
from Ui_NoLoadSimulation import  Ui_MainWindow
from PyQt5 import QtCore, QtGui, QtWidgets

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow) 
    MainWindow.show()
    sys.exit(app.exec_())