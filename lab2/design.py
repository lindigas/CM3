# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'design.ui'
#
# Created by: PyQt5 UI code generator 5.11.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(742, 585)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setContentsMargins(5, 5, 5, 5)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.tableWidget = QtWidgets.QTableWidget(self.centralwidget)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.gridLayout_2.addWidget(self.tableWidget, 1, 0, 1, 2)
        self.textEdit = QtWidgets.QTextEdit(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.textEdit.setFont(font)
        self.textEdit.setObjectName("textEdit")
        self.gridLayout_2.addWidget(self.textEdit, 0, 1, 1, 1)
        self.gridLayout_4 = QtWidgets.QGridLayout()
        self.gridLayout_4.setSizeConstraint(QtWidgets.QLayout.SetMinimumSize)
        self.gridLayout_4.setContentsMargins(5, 5, 5, 5)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setObjectName("label_5")
        self.gridLayout_4.addWidget(self.label_5, 0, 0, 1, 1)
        self.x_num_2 = QtWidgets.QLineEdit(self.centralwidget)
        self.x_num_2.setObjectName("x_num_2")
        self.gridLayout_4.addWidget(self.x_num_2, 0, 1, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setObjectName("label_6")
        self.gridLayout_4.addWidget(self.label_6, 1, 0, 1, 1)
        self.y_num_2 = QtWidgets.QLineEdit(self.centralwidget)
        self.y_num_2.setObjectName("y_num_2")
        self.gridLayout_4.addWidget(self.y_num_2, 1, 1, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setObjectName("label_7")
        self.gridLayout_4.addWidget(self.label_7, 2, 0, 1, 1)
        self.accuracy_2 = QtWidgets.QLineEdit(self.centralwidget)
        self.accuracy_2.setObjectName("accuracy_2")
        self.gridLayout_4.addWidget(self.accuracy_2, 2, 1, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setObjectName("label_8")
        self.gridLayout_4.addWidget(self.label_8, 3, 0, 1, 1)
        self.breaks_2 = QtWidgets.QLineEdit(self.centralwidget)
        self.breaks_2.setObjectName("breaks_2")
        self.gridLayout_4.addWidget(self.breaks_2, 3, 1, 1, 1)
        self.btnTest_2 = QtWidgets.QPushButton(self.centralwidget)
        self.btnTest_2.setObjectName("btnTest_2")
        self.gridLayout_4.addWidget(self.btnTest_2, 4, 0, 1, 1)
        self.btnMain_2 = QtWidgets.QPushButton(self.centralwidget)
        self.btnMain_2.setObjectName("btnMain_2")
        self.gridLayout_4.addWidget(self.btnMain_2, 4, 1, 1, 1)
        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_2.setObjectName("pushButton_2")
        self.gridLayout_4.addWidget(self.pushButton_2, 5, 0, 1, 1)
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout_4.addWidget(self.pushButton, 5, 1, 1, 1)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.gridLayout_4.addLayout(self.gridLayout, 6, 1, 1, 1)
        self.gridLayout_2.addLayout(self.gridLayout_4, 0, 0, 1, 1)
        self.gridLayout_3.addLayout(self.gridLayout_2, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 742, 26))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label_5.setText(_translate("MainWindow", "Разбиений по х"))
        self.x_num_2.setText(_translate("MainWindow", "16"))
        self.label_6.setText(_translate("MainWindow", "Разбиений по у"))
        self.y_num_2.setText(_translate("MainWindow", "16"))
        self.label_7.setText(_translate("MainWindow", "Точность метода"))
        self.accuracy_2.setText(_translate("MainWindow", "0.00000001"))
        self.label_8.setText(_translate("MainWindow", "Ограничение шагов"))
        self.breaks_2.setText(_translate("MainWindow", "1000"))
        self.btnTest_2.setText(_translate("MainWindow", "Решить тестовую задачу"))
        self.btnMain_2.setText(_translate("MainWindow", "Решить основную задачу"))
        self.pushButton_2.setText(_translate("MainWindow", "Показать тестовую задачу"))
        self.pushButton.setText(_translate("MainWindow", "Показать основную задачу"))

