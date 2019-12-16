import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QWidget, QTableWidget, QTableWidgetItem, QHeaderView
import numpy as np
import design
import test
import maintask
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import math

class MainWindow(QtWidgets.QMainWindow, design.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.btnTest_2.clicked.connect(self.TestTask)
        self.btnMain_2.clicked.connect(self.MainTask)
        self.pushButton_2.clicked.connect(self.showTest)
        self.pushButton.clicked.connect(self.showMain)

    def showMain(self):
        plt.imshow(mpimg.imread('main.png'))
        plt.show()

    def showTest(self):
        plt.imshow(mpimg.imread('test.png'))
        plt.show()
    
    def TestTask(self):
        text_list = ["При решении Р.С. с помощью метода Зейделя \n Nmax = ", " и eps = ", 
                     " за S = ", " итераций было получено решение \n с точностью eps_max = ", 
                     " и максимальной невязкой ||r|| = ","\nМаксимальная погрешность решения тестовой задачи = ", " при х = ",
                                " и у = "]
        n = np.int(self.x_num_2.text())
        m = np.int(self.y_num_2.text())
        h = 2./n
        k = 1./m
        eps = np.float64(self.accuracy_2.text())
        nmax = np.int(self.breaks_2.text())
        x, y, v = test.grid(0, 2, 0, 1, n, m)
        ss, ee, v = test.zeidel(n, m, 0, 2, 0, 1, x, y, v, eps, nmax)
        max_nev = test.discrepancy(v,n,m,x,y,h,k)
        self.Table(v, n, m)
        Max, MaxX, MaxY = test.accuracy(n, m, x, y, v)
        text_list.insert(1,str(nmax))
        text_list.insert(3,str(eps))
        text_list.insert(5,str(ss))
        text_list.insert(7,str(ee))
        text_list.insert(9,str(max_nev))
        text_list.insert(11,str(Max))
        text_list.insert(13,str(MaxX))
        text_list.insert(15,str(MaxY))
        text = ''.join(text_list)
        self.textEdit.setText(text)


    def MainTask(self):
        text_list = ["При решении Р.С. с помощью метода Зейделя Nmax = ", ", eps = ", 
                     "\n За S = ", " итераций было получено решение с точностью eps_max = ", 
                     " и максимальной невязкой ||r|| = "]
        tail = ["Максимальная погрешность решения основной задачи = ", " при х = ",
                                " и у = "]
        n = np.int(self.x_num_2.text())
        m = np.int(self.y_num_2.text())
        h = 2./n
        k = 1./m
        eps = np.float64(self.accuracy_2.text())
        nmax = np.int(self.breaks_2.text())
        x1, y1, v1, x2, y2, v2 = maintask.gridd(0, 2, 0, 1, n, m)
        ss, ee, v1 = maintask.zeidel(n, m, 0, 2, 0, 1, x1, y1, v1, eps, nmax)
        max_nev = maintask.discrepancy(v1, n, m, x1, y1, h, k)
        ss2, ee2, v2 = maintask.zeidel(2 * n, 2 * m, 0, 2, 0, 1, x2, y2, v2, eps, nmax)
        max_nev2 = maintask.discrepancy(v2, 2*n, 2*m, x2, y2, h/2., k/2.)
        self.Table(v1, n, m)
        Max, MaxX, MaxY = maintask.accuracy(n, m, x1, y1, v1, v2)
        text_list.insert(1,str(nmax))
        text_list.insert(3,str(eps))
        text_list.insert(5,str(ss))
        text_list.insert(7,str(ee))
        text_list.insert(9,str(max_nev))
        text = ''.join(text_list)
        text_list[1] = str(nmax * 2)
        text_list[3] = str(eps)
        text_list[5] = str(ss2) 
        text_list[7] = str(ee2)
        text_list[9] = str(max_nev2)
        text_list.insert(0,"\n \n")
        text1 = ''.join(text_list)
        tail.insert(1,str(Max))
        tail.insert(3,str(MaxX))
        tail.insert(5,str(MaxY))
        tail.insert(0,"\n \n")
        text2 = ''.join(tail)
        self.textEdit.setText(text+text1+text2)


    def Table(self, v, n, m):
        coef = 1
        if ((n>50) and (m>50)):
            coef = 5
        if ((n>100) and (m>100)):
            coef = 10
        if (n>300) and (m>300):
            coef = 30
        self.tableWidget.setRowCount(0)
        self.tableWidget.setRowCount(int(m/coef))
        self.tableWidget.setColumnCount(int(n/coef))
        for i in range(0,int(n/coef)):
            for j in range(0, int(m/coef)):
                ii = math.floor(i*coef)
                jj = math.floor(j*coef)
                self.tableWidget.setItem(i, j, QTableWidgetItem(str(v[ii][jj])))
                header = self.tableWidget.horizontalHeader()
                for k in range(0, jj):
                    header.setSectionResizeMode(k, QtWidgets.QHeaderView.ResizeToContents)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)

    main = MainWindow()
    main.show()

    sys.exit(app.exec_())
