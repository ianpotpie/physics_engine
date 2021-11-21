from PyQt5 import QtGui
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtGui import QPainter, QBrush, QPen
from PyQt5.QtCore import Qt
import sys
from PyQt5 import QtWidgets, QtCore

#
# class Window(QMainWindow):
#     def __init__(self):
#         super().__init__()
#
#         self.title = "PyQt5 Drawing Rectangle"
#         self.top = 100
#         self.left = 100
#         self.width = 680
#         self.height = 500
#
#         self.InitWindow()
#
#     def InitWindow(self):
#         self.setWindowIcon(QtGui.QIcon("icon.png"))
#         self.setWindowTitle(self.title)
#         self.setGeometry(self.top, self.left, self.width, self.height)
#         self.show()
#
#     def paintEvent(self, e):
#         painter = QPainter(self)
#         painter.setPen(QPen(Qt.black, 5, Qt.SolidLine))
#         # painter.setBrush(QBrush(Qt.red, Qt.SolidPattern))
#         # painter.setBrush(QBrush(Qt.green))
#         painter.drawRect(100, 15, 400, 200)

class MyWidget(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setGeometry(30,30,600,400)
        self.begin = QtCore.QPoint()
        self.end = QtCore.QPoint()
        self.show()

    def paintEvent(self, event):
        qp = QtGui.QPainter(self)
        br = QtGui.QBrush(QtGui.QColor(100, 10, 10, 40))
        qp.setBrush(br)
        qp.drawRect(QtCore.QRect(self.begin, self.end))

    def mousePressEvent(self, event):
        self.begin = event.pos()
        self.end = event.pos()
        self.update()

    def mouseMoveEvent(self, event):
        self.end = event.pos()
        self.update()

    def mouseReleaseEvent(self, event):
        self.begin = event.pos()
        self.end = event.pos()
        self.update()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = MyWidget()
    window.show()
    app.aboutToQuit.connect(app.deleteLater)
    sys.exit(app.exec_())

# App = QApplication(sys.argv)
# window = Window()
# sys.exit(App.exec())
