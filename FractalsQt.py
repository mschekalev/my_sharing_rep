from PySide2 import QtCore, QtGui, QtWidgets
import sys
import math


class MyWidget(QtWidgets.QWidget):
    def Coch(self, pnt, n, x1, y1, x2, y2):
        if n == 0:
            pnt.drawLine(x1, self.size().height() - y1, x2, self.size().height() - y2)
        else:
            f = math.pi / 3
            dx, dy = (x2 - x1) / 3, (y2 - y1) / 3
            ax, ay = x1 + dx, y1 + dy
            bx, by = x2 - dx, y2 - dy
            cx, cy = ax + dx * math.cos(f) - dy * math.sin(f), ay + dx * math.sin(f) + dy * math.cos(f)
            self.Coch(pnt, n - 1, x1, y1, ax, ay)
            self.Coch(pnt, n - 1, ax, ay, cx, cy)
            self.Coch(pnt, n - 1, cx, cy, bx, by)
            self.Coch(pnt, n - 1, bx, by, x2, y2)
        
    def CochHex(self, pnt, n, a, x, y):
        self.Coch(pnt, n, x, y, x, y + a)
        self.Coch(pnt, n, x, y + a, x + math.sqrt(3) * a / 2, y + a / 2)
        self.Coch(pnt, n, x + math.sqrt(3) * a / 2, y + a / 2, x, y)
    
    def Triangle(self, pnt, n, a, x, y):
        if n == 0:
            pnt.drawLine(x, y, x + a, y)
            pnt.drawLine(x, y, x + a / 2, y - math.sqrt(3) * a / 2)
            pnt.drawLine(x + a / 2, y - math.sqrt(3) * a / 2, x + a, y)
        else:
            self.Triangle(pnt, n - 1, a / 2, x, y)
            self.Triangle(pnt, n - 1, a / 2, x + a / 4, y - math.sqrt(3) * a / 4)
            self.Triangle(pnt, n - 1, a / 2, x + a / 2, y)
            
    def Minkowski(self, pnt, n, xa, ya, xi, yi):
        if n == 0:            
            pnt.drawLine(xa, ya, xi, yi)
        else:
            dx, dy = xi - xa, yi - ya
            xb, yb = xa + dx / 4, ya + dy / 4
            xe, ye = xa + dx / 2, ya + dy / 2
            xh, yh = xa + 3 * dx / 4, ya + 3 * dy / 4
            xc, yc = xb + ye - yb, yb - xe + xb
            xd, yd = xc + xe - xb, yc + ye - yb
            xf, yf = xe - yh + ye, ye + xh - xe
            xg, yg = xf + xh - xe, yf + yh - ye
            
            self.Minkowski(pnt, n - 1, xa, ya, xb, yb)
            self.Minkowski(pnt, n - 1, xb, yb, xc, yc)
            self.Minkowski(pnt, n - 1, xc, yc, xd, yd)
            self.Minkowski(pnt, n - 1, xd, yd, xe, ye)
            self.Minkowski(pnt, n - 1, xe, ye, xf, yf)
            self.Minkowski(pnt, n - 1, xf, yf, xg, yg)
            self.Minkowski(pnt, n - 1, xg, yg, xh, yh)
            self.Minkowski(pnt, n - 1, xh, yh, xi, yi)            
     
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.__n, self.__type = 0, "Coch"
        self.resize(300, 300)
        self.show()
    
    def paintEvent(self, event):
        painter = QtGui.QPainter()
        painter.begin(self)
        X, Y = self.size().width(), self.size().height()
        S = 0
        if self.__type == "Coch":
            a = 3 * X / (2 * math.sqrt(3))
            if a > Y - 1:
                a = Y - 1
            self.CochHex(painter, self.__n, a, (X - math.sqrt(3) * a / 3) / 2, (Y - a) / 2)
        elif self.__type == "Minkowski line":
            if X / 3 >= Y / 2:
                S = 1.5 * (X / 3 - Y / 2)
            self.Minkowski(painter, self.__n, S, Y / 2, X - S, Y / 2)
        else:
            h = math.sqrt(3) * X / 2
            if h >= Y - 1:
                S = 1.15 * (h - Y)
            self.Triangle(painter, self.__n, X - S, S / 2, Y / 2 + math.sqrt(3) * (X - S) / 4)
        painter.end()
        
    def setValue(self, val):
        self.__n = val
        self.repaint()
    
    def setFractalType(self, val):
        self.__type = val
        self.repaint()


class MyWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self.resize(300, 300)
        self.setMinimumSize(300, 300)
        self.setWindowTitle("Painter demo")
        self.Widget = MyWidget(self)
        self.Widget.setGeometry(0, 30, 300, 300)
        
    def resizeEvent(self, event):
        self.Widget.setGeometry(0, 0, self.width(), self.height())


app = QtWidgets.QApplication(sys.argv)
Window = MyWindow()
spin = QtWidgets.QSpinBox(Window)
spin.setGeometry(10, 10, 100, 30)
combo = QtWidgets.QComboBox(Window)
combo.setGeometry(130, 10, 130, 30)
combo.addItems('Coch,Minkowski line,Sierpinski triangle'.split(','))
QtCore.QObject.connect(spin, QtCore.SIGNAL("valueChanged(int)"),
                       Window.Widget.setValue)
QtCore.QObject.connect(combo, QtCore.SIGNAL("currentIndexChanged(QString)"),
                       Window.Widget.setFractalType)
Window.show()
app.exec_()