from PySide2 import QtCore, QtWidgets
import sys
import random


class Game(QtCore.QObject):
    def __init__(self):
        QtCore.QObject.__init__(self)
        self.__initialization()
        self.__ends = 0
        
    def __initialization(self):
        self.__num = randomNum()
        print(self.__num)
        self.__usnum, self.__turns, self.__prev = '0', 0, '0'
        
    def endgame(self):
        self.__ends += 1
        if self.__ends % 2 == 1:
            self.emit(QtCore.SIGNAL("changedText(QString)"), '<center><b><font size="7"> You surrendered.</font></b><center>')
            self.emit(QtCore.SIGNAL("changedComment(QString)"), '<center><i><font size="4"> Answer: <b>%s</b> </font></i><center>' % str(self.__num))
            self.emit(QtCore.SIGNAL("endedGame(QString)"), 'New Game')
            self.emit(QtCore.SIGNAL("endedGame(bool)"), False)
        else:
            self.emit(QtCore.SIGNAL("changedText(QString)"), '')
            self.emit(QtCore.SIGNAL("changedComment(QString)"), '')
            self.emit(QtCore.SIGNAL("endedGame(QString)"), 'Give up')
            self.emit(QtCore.SIGNAL("changedLine(QString)"), '')
            self.emit(QtCore.SIGNAL("endedGame(bool)"), True)
            self.__initialization()
        
    def setNum(self, value):
        self.__usnum = value
               
    def check(self):
        num = self.__usnum
        if not num.isdigit() or int(num) < 1023 or num[0] == '0' or len(set(num)) != 4 or len(num) > 4:
            self.emit(QtCore.SIGNAL("changedText(QString)"), '<center><b><font size="7" color="red">Wrong value!</font></b><center>')
            self.emit(QtCore.SIGNAL("changedLine(QString)"), '')
        else:
            if num != self.__prev:
                self.__turns += 1
                self.__prev = num
            cows, bulls = 0, 0
            for i in range(4):
                if num[i] in self.__num:
                    if i == self.__num.find(num[i]):
                        bulls += 1
                    else:
                        cows += 1
            if bulls == 4:
                self.__ends += 1
                self.emit(QtCore.SIGNAL("changedText(QString)"), '<center><b><font size="7" color="blue">You won!</font></i><center>')
                self.emit(QtCore.SIGNAL("changedComment(QString)"), '<center><i><font size="4"> Turns made: <b>%s</b></font></i><center>' % str(self.__turns))
                self.emit(QtCore.SIGNAL("endedGame(QString)"), 'New Game')
                self.emit(QtCore.SIGNAL("endedGame(bool)"), False)
            else:
                self.emit(QtCore.SIGNAL("changedText(QString)"), '<center><b><font size="7">Bulls: %s; Cows: %s</font></b><center>' % (str(bulls), str(cows)))


def randomNum():
    res = ''
    List = [i for i in range(1, 10)]
    curr = random.choice(List)
    res += str(curr)
    List.pop(curr - 1)
    List = [0] + List
    for i in range(3):
        curr = random.choice(List)
        res += str(curr)
        List.pop(List.index(curr))
    return res
        
  
app = QtWidgets.QApplication(sys.argv)
window = QtWidgets.QMainWindow()
window.setWindowTitle('Bulls n Cows')
window.setMinimumSize(265, 130)
window.setMaximumSize(265, 130)
enterLable = QtWidgets.QLabel('<b><font size="4" >Enter number:</font></b>', window)
enterLable.setGeometry(10, 10, 100, 30)
LineEdit = QtWidgets.QLineEdit(window)
LineEdit.setGeometry(115, 10, 90, 30)
Button = QtWidgets.QPushButton('Check', window)
Button.setGeometry(210, 9, 50, 31)
resLable = QtWidgets.QLabel(window)
resLable.setGeometry(10, 50, 255, 40)
comLable = QtWidgets.QLabel(window)
comLable.setGeometry(90, 95, 175, 35)
sButton = QtWidgets.QPushButton('Give up', window)
sButton.setGeometry(10, 90, 70, 35)

MyGame = Game()
QtCore.QObject.connect(LineEdit, QtCore.SIGNAL("textChanged(QString)"),
                       MyGame.setNum)
QtCore.QObject.connect(Button, QtCore.SIGNAL("clicked()"),
                       MyGame.check)
QtCore.QObject.connect(sButton, QtCore.SIGNAL("clicked()"),
                       MyGame.endgame)
QtCore.QObject.connect(MyGame, QtCore.SIGNAL("changedText(QString)"), 
                       resLable.setText)
QtCore.QObject.connect(MyGame, QtCore.SIGNAL("changedComment(QString)"), 
                       comLable.setText)
QtCore.QObject.connect(MyGame, QtCore.SIGNAL("endedGame(QString)"), 
                       sButton.setText)
QtCore.QObject.connect(MyGame, QtCore.SIGNAL("endedGame(bool)"), 
                       Button.setEnabled)
QtCore.QObject.connect(MyGame, QtCore.SIGNAL("changedLine(QString)"), 
                       LineEdit.setText)
window.show()
app.exec_()
