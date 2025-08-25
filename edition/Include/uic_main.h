/********************************************************************************
** Form generated from reading UI file 'main.ui'
**
** Created by: Qt User Interface Compiler version 6.2.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UIC_MAIN_H
#define UIC_MAIN_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtOpenGLWidgets/QOpenGLWidget>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStackedWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Assets
{
public:
    QAction *actionExit;
    QAction *actionExport_obj;
    QAction *actionExport_ma;
    QAction *actionLoad_ma;
    QAction *actionGL;
    QAction *actionVectorized;
    QAction *actionExport_svg;
    QAction *actionExport_mi;
    QAction *action128x128;
    QAction *action256x256;
    QAction *action512x512;
    QAction *action1024x1024;
    QAction *action2048x2048;
    QWidget *centralwidget;
    QHBoxLayout *hboxLayout;
    QWidget *Objects_widget;
    QVBoxLayout *verticalLayout;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_2;
    QPushButton *btn_addDetails;
    QPushButton *btn_loadShader;
    QPushButton *openGaussianFile;
    QPushButton *btn_saveGaussian;
    QPushButton *btn_exporthighres;
    QPushButton *btn_openTerrain;
    QPushButton *btn_clear;
    QSpacerItem *verticalSpacer;
    QComboBox *combo_exporthighres;
    QGroupBox *groupBox_3;
    QCheckBox *check_showInfluence;
    QSlider *slider_nbgaussians;
    QLabel *label_nbgaussians;
    QLabel *label_3;
    QLabel *label_9;
    QLabel *label_7;
    QSlider *slider_amplitude;
    QSlider *slider_noiseLevel;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout_3;
    QPushButton *btn_graphTool;
    QPushButton *btn_handTool;
    QPushButton *btn_eraseTool;
    QPushButton *btn_moveTool;
    QStackedWidget *toolOptionsWidget;
    QWidget *page_2;
    QWidget *page_3;
    QGroupBox *groupBox_brush;
    QSlider *slider_brushthreshold;
    QLabel *label_4;
    QPushButton *btn_applyTool;
    QGroupBox *groupBox_4;
    QListWidget *list_brushes;
    QPushButton *btn_openBrush;
    QPushButton *btn_saveBrush;
    QWidget *page_4;
    QGroupBox *groupBox_graph;
    QSlider *slider_depthGraph;
    QLabel *label_5;
    QLabel *label_8;
    QSlider *slider_blendThresholdGraph;
    QCheckBox *check_influence;
    QWidget *page;
    QGroupBox *groupBox_edit;
    QRadioButton *radio_erase;
    QRadioButton *radio_amplitude;
    QRadioButton *radio_warp;
    QRadioButton *radio_move;
    QRadioButton *radio_amplitude_lr;
    QRadioButton *radio_amplitude_hr;
    QCheckBox *check_saveLogs;
    QPushButton *btn_curveTool;
    QPushButton *btn_AdrienGraph;
    QPushButton *btn_test;
    QOpenGLWidget *raytracingwidget;
    QMenuBar *menuBar;
    QMenu *menuOptions;
    QMenu *menuRender_resolution;

    void setupUi(QMainWindow *Assets)
    {
        if (Assets->objectName().isEmpty())
            Assets->setObjectName(QString::fromUtf8("Assets"));
        Assets->resize(1462, 964);
        Assets->setMinimumSize(QSize(420, 300));
        Assets->setStyleSheet(QString::fromUtf8("/*-----QWidget-----*/\n"
"QWidget\n"
"{\n"
"	background-color: #3a3a3a;\n"
"	color: #fff;\n"
"	selection-background-color: #b78620;\n"
"	selection-color: #000;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QLabel-----*/\n"
"QLabel\n"
"{\n"
"	background-color: transparent;\n"
"	color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QMenuBar-----*/\n"
"QMenuBar \n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(57, 57, 57, 255),stop:1 rgba(50, 50, 50, 255));\n"
"	border: 1px solid #000;\n"
"	color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"QMenuBar::item \n"
"{\n"
"	background-color: transparent;\n"
"\n"
"}\n"
"\n"
"\n"
"QMenuBar::item:selected \n"
"{\n"
"	background-color: rgba(183, 134, 32, 20%);\n"
"	border: 1px solid #b78620;\n"
"	color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"QMenuBar::item:pressed \n"
"{\n"
"	background-color: rgb(183, 134, 32);\n"
"	border: 1px solid #b78620;\n"
"	color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QMenu-----*/\n"
"QMenu\n"
"{\n"
"    background-color: qlineargradient(spread:r"
                        "epeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(57, 57, 57, 255),stop:1 rgba(50, 50, 50, 255));\n"
"    border: 1px solid #222;\n"
"    padding: 4px;\n"
"	color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"QMenu::item\n"
"{\n"
"    background-color: transparent;\n"
"    padding: 2px 20px 2px 20px;\n"
"\n"
"}\n"
"\n"
"\n"
"QMenu::separator\n"
"{\n"
"   	background-color: rgb(183, 134, 32);\n"
"	height: 1px;\n"
"\n"
"}\n"
"\n"
"\n"
"QMenu::item:disabled\n"
"{\n"
"    color: #555;\n"
"    background-color: transparent;\n"
"    padding: 2px 20px 2px 20px;\n"
"\n"
"}\n"
"\n"
"\n"
"QMenu::item:selected\n"
"{\n"
"	background-color: rgba(183, 134, 32, 20%);\n"
"	border: 1px solid #b78620;\n"
"	color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QToolBar-----*/\n"
"QToolBar\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(69, 69, 69, 255),stop:1 rgba(58, 58, 58, 255));\n"
"	border-top: none;\n"
"	border-bottom: 1px solid #4f4f4f;\n"
"	border-left: 1px solid #4f4f4f;\n"
"	border-right: 1px solid #4"
                        "f4f4f;\n"
"\n"
"}\n"
"\n"
"\n"
"QToolBar::separator\n"
"{\n"
"	background-color: #2e2e2e;\n"
"	width: 1px;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QToolButton-----*/\n"
"QToolButton \n"
"{\n"
"	background-color: transparent;\n"
"	color: #fff;\n"
"	padding: 5px;\n"
"	padding-left: 8px;\n"
"	padding-right: 8px;\n"
"	margin-left: 1px;\n"
"}\n"
"\n"
"\n"
"QToolButton:hover\n"
"{\n"
"	background-color: rgba(183, 134, 32, 20%);\n"
"	border: 1px solid #b78620;\n"
"	color: #fff;\n"
"	\n"
"}\n"
"\n"
"\n"
"QToolButton:pressed\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(57, 57, 57, 255),stop:1 rgba(50, 50, 50, 255));\n"
"	border: 1px solid #b78620;\n"
"\n"
"}\n"
"\n"
"\n"
"QToolButton:checked\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(57, 57, 57, 255),stop:1 rgba(50, 50, 50, 255));\n"
"	border: 1px solid #222;\n"
"}\n"
"\n"
"\n"
"/*-----QPushButton-----*/\n"
"QPushButton\n"
"{\n"
"	background-color: qlineargradient(spread:rep"
                        "eat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(84, 84, 84, 255),stop:1 rgba(59, 59, 59, 255));\n"
"	color: #ffffff;\n"
"	min-width: 80px;\n"
"	border-style: solid;\n"
"	border-width: 1px;\n"
"	border-radius: 3px;\n"
"	border-color: #051a39;\n"
"	padding: 5px;\n"
"\n"
"}\n"
"\n"
"\n"
"QPushButton::flat\n"
"{\n"
"	background-color: transparent;\n"
"	border: none;\n"
"	color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"QPushButton::disabled\n"
"{\n"
"	background-color: #404040;\n"
"	color: #656565;\n"
"	border-color: #051a39;\n"
"\n"
"}\n"
"\n"
"\n"
"QPushButton::hover\n"
"{\n"
"	background-color: rgba(183, 134, 32, 20%);\n"
"	border: 1px solid #b78620;\n"
"\n"
"}\n"
"\n"
"\n"
"QPushButton::pressed\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(74, 74, 74, 255),stop:1 rgba(49, 49, 49, 255));\n"
"	border: 1px solid #b78620;\n"
"\n"
"}\n"
"\n"
"\n"
"QPushButton::checked\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(74, 74, 74, 255),stop"
                        ":1 rgba(49, 49, 49, 255));\n"
"	border: 1px solid #222;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QLineEdit-----*/\n"
"QLineEdit\n"
"{\n"
"	background-color: #131313;\n"
"	color : #eee;\n"
"	border: 1px solid #343434;\n"
"	border-radius: 2px;\n"
"	padding: 3px;\n"
"	padding-left: 5px;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QPlainTExtEdit-----*/\n"
"QPlainTextEdit\n"
"{\n"
"	background-color: #131313;\n"
"	color : #eee;\n"
"	border: 1px solid #343434;\n"
"	border-radius: 2px;\n"
"	padding: 3px;\n"
"	padding-left: 5px;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QPlainTExtEdit-----*/\n"
"QTextEdit\n"
"{\n"
"	background-color: #131313;\n"
"	color : #eee;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QTabBar-----*/\n"
"QTabBar::tab\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(84, 84, 84, 255),stop:1 rgba(59, 59, 59, 255));\n"
"	color: #ffffff;\n"
"	border-style: solid;\n"
"	border-width: 1px;\n"
"	border-color: #666;\n"
"	border-bottom: none;\n"
"	padding: 5px;\n"
"	padding-left: 15px;\n"
"	padding-right: "
                        "15px;\n"
"\n"
"}\n"
"\n"
"\n"
"QTabWidget::pane \n"
"{\n"
"	background-color: red;\n"
"	border: 1px solid #666;\n"
"	top: 1px;\n"
"\n"
"}\n"
"\n"
"\n"
"QTabBar::tab:last\n"
"{\n"
"	margin-right: 0; \n"
"\n"
"}\n"
"\n"
"\n"
"QTabBar::tab:first:!selected\n"
"{\n"
"	background-color: #0c0c0d;\n"
"	margin-left: 0px;\n"
"\n"
"}\n"
"\n"
"\n"
"QTabBar::tab:!selected\n"
"{\n"
"	color: #b1b1b1;\n"
"	border-bottom-style: solid;\n"
"	background-color: #0c0c0d;\n"
"\n"
"}\n"
"\n"
"\n"
"QTabBar::tab:selected\n"
"{\n"
"	margin-bottom: 0px;\n"
"\n"
"}\n"
"\n"
"\n"
"QTabBar::tab:!selected:hover\n"
"{\n"
"	border-top-color: #b78620;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QComboBox-----*/\n"
"QComboBox\n"
"{\n"
"    background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(84, 84, 84, 255),stop:1 rgba(59, 59, 59, 255));\n"
"    border: 1px solid #000;\n"
"    padding-left: 6px;\n"
"    color: #ffffff;\n"
"    height: 20px;\n"
"\n"
"}\n"
"\n"
"\n"
"QComboBox::disabled\n"
"{\n"
"	background-color: #404040;\n"
""
                        "	color: #656565;\n"
"	border-color: #051a39;\n"
"\n"
"}\n"
"\n"
"\n"
"QComboBox:on\n"
"{\n"
"    background-color: #b78620;\n"
"	color: #000;\n"
"\n"
"}\n"
"\n"
"\n"
"QComboBox QAbstractItemView\n"
"{\n"
"    background-color: #383838;\n"
"    color: #ffffff;\n"
"    border: 1px solid black;\n"
"    selection-background-color: #b78620;\n"
"    outline: 0;\n"
"\n"
"}\n"
"\n"
"\n"
"QComboBox::drop-down\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(57, 57, 57, 255),stop:1 rgba(50, 50, 50, 255));\n"
"    subcontrol-origin: padding;\n"
"    subcontrol-position: top right;\n"
"    width: 15px;\n"
"    border-left-width: 1px;\n"
"    border-left-color: black;\n"
"    border-left-style: solid; \n"
"\n"
"}\n"
"\n"
"\n"
"QComboBox::down-arrow\n"
"{\n"
"    image: url(://images/arrow-down.png);\n"
"    width: 8px;\n"
"    height: 8px;\n"
"}\n"
"\n"
"\n"
"/*-----QSpinBox & QDateTimeEdit-----*/\n"
"QSpinBox,\n"
"QDoubleSpinBox,\n"
"QDateTimeEdit \n"
"{\n"
"    background-co"
                        "lor: #131313;\n"
"	color : #eee;\n"
"	border: 1px solid #343434;\n"
"	padding: 3px;\n"
"	padding-left: 5px;\n"
"    border-radius : 2px;\n"
"\n"
"}\n"
"\n"
"\n"
"QSpinBox::up-button, \n"
"QDoubleSpinBox::up-button,\n"
"QDateTimeEdit::up-button\n"
"{\n"
"	border-top-right-radius:2px;\n"
"	background-color: #777777;\n"
"    width: 16px; \n"
"    border-width: 1px;\n"
"\n"
"}\n"
"\n"
"\n"
"QSpinBox::up-button:hover, \n"
"QDoubleSpinBox:hover,\n"
"QDateTimeEdit::up-button:hover\n"
"{\n"
"	background-color: #585858;\n"
"\n"
"}\n"
"\n"
"\n"
"QSpinBox::up-button:pressed, \n"
"QDoubleSpinBox::up-button:pressed, \n"
"QDateTimeEdit::up-button:pressed\n"
"{\n"
"	background-color: #252525;\n"
"    width: 16px; \n"
"    border-width: 1px;\n"
"\n"
"}\n"
"\n"
"\n"
"QSpinBox::up-arrow,\n"
"QDoubleSpinBox::up-arrow,\n"
"QDateTimeEdit::up-arrow\n"
"{\n"
"    image: url(:/images/arrow-up.png);\n"
"    width: 7px;\n"
"    height: 7px;\n"
"\n"
"}\n"
"\n"
"\n"
"QSpinBox::down-button, \n"
"QDoubleSpinBox::down-button, \n"
"QDateTime"
                        "Edit::down-button\n"
"{\n"
"	border-bottom-right-radius:2px;\n"
"	background-color: #777777;\n"
"    width: 16px; \n"
"    border-width: 1px;\n"
"\n"
"}\n"
"\n"
"\n"
"QSpinBox::down-button:hover, \n"
"QDoubleSpinBox::down-button:hover, \n"
"QDateTimeEdit::down-button:hover\n"
"{\n"
"	background-color: #585858;\n"
"\n"
"}\n"
"\n"
"\n"
"QSpinBox::down-button:pressed, \n"
"QDoubleSpinBox::down-button:pressed, \n"
"QDateTimeEdit::down-button:pressed\n"
"{\n"
"	background-color: #252525;\n"
"    width: 16px; \n"
"    border-width: 1px;\n"
"\n"
"}\n"
"\n"
"\n"
"QSpinBox::down-arrow,\n"
"QDoubleSpinBox::down-arrow,\n"
"QDateTimeEdit::down-arrow\n"
"{\n"
"    image: url(:/images/arrow-down.png);\n"
"    width: 7px;\n"
"    height: 7px;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QGroupBox-----*/\n"
"QGroupBox \n"
"{\n"
"    border: 1px solid;\n"
"    border-color: #666666;\n"
"	border-radius: 5px;\n"
"    margin-top: 25px;\n"
"\n"
"}\n"
"\n"
"\n"
"QGroupBox::title  \n"
"{\n"
"    background-color: transparent;\n"
"    color: #eee"
                        ";\n"
"    subcontrol-origin: margin;\n"
"    padding: 5px;\n"
"	border-top-left-radius: 3px;\n"
"	border-top-right-radius: 3px;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QHeaderView-----*/\n"
"QHeaderView::section\n"
"{\n"
"    background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(60, 60, 60, 255),stop:1 rgba(50, 50, 50, 255));\n"
"	border: 1px solid #000;\n"
"    color: #fff;\n"
"    text-align: left;\n"
"	padding: 4px;\n"
"	\n"
"}\n"
"\n"
"\n"
"QHeaderView::section:disabled\n"
"{\n"
"    background-color: #525251;\n"
"    color: #656565;\n"
"\n"
"}\n"
"\n"
"\n"
"QHeaderView::section:checked\n"
"{\n"
"    background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(60, 60, 60, 255),stop:1 rgba(50, 50, 50, 255));\n"
"    color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"QHeaderView::section::vertical::first,\n"
"QHeaderView::section::vertical::only-one\n"
"{\n"
"    border-top: 1px solid #353635;\n"
"\n"
"}\n"
"\n"
"\n"
"QHeaderView::section::vertical\n"
"{\n"
"    border-top: "
                        "1px solid #353635;\n"
"\n"
"}\n"
"\n"
"\n"
"QHeaderView::section::horizontal::first,\n"
"QHeaderView::section::horizontal::only-one\n"
"{\n"
"    border-left: 1px solid #353635;\n"
"\n"
"}\n"
"\n"
"\n"
"QHeaderView::section::horizontal\n"
"{\n"
"    border-left: 1px solid #353635;\n"
"\n"
"}\n"
"\n"
"\n"
"QTableCornerButton::section\n"
"{\n"
"    background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(60, 60, 60, 255),stop:1 rgba(50, 50, 50, 255));\n"
"	border: 1px solid #000;\n"
"    color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QTreeWidget-----*/\n"
"QTreeView\n"
"{\n"
"	show-decoration-selected: 1;\n"
"	alternate-background-color: #3a3a3a;\n"
"	selection-color: #fff;\n"
"	background-color: #2d2d2d;\n"
"	border: 1px solid gray;\n"
"	padding-top : 5px;\n"
"	color: #fff;\n"
"	font: 8pt;\n"
"\n"
"}\n"
"\n"
"\n"
"QTreeView::item:selected\n"
"{\n"
"	color:#fff;\n"
"	background-color: #b78620;\n"
"	border-radius: 0px;\n"
"\n"
"}\n"
"\n"
"\n"
"QTreeView::item:!selected:hover\n"
"{\n"
"  "
                        "  background-color: #262626;\n"
"    border: none;\n"
"    color: white;\n"
"\n"
"}\n"
"\n"
"\n"
"QTreeView::branch:has-children:!has-siblings:closed,\n"
"QTreeView::branch:closed:has-children:has-siblings \n"
"{\n"
"	image: url(://tree-closed.png);\n"
"\n"
"}\n"
"\n"
"\n"
"QTreeView::branch:open:has-children:!has-siblings,\n"
"QTreeView::branch:open:has-children:has-siblings  \n"
"{\n"
"	image: url(://tree-open.png);\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QListView-----*/\n"
"QListView \n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(83, 83, 83, 255),stop:0.293269 rgba(81, 81, 81, 255),stop:0.634615 rgba(79, 79, 79, 255),stop:1 rgba(83, 83, 83, 255));\n"
"    border : none;\n"
"    color: white;\n"
"    show-decoration-selected: 1; \n"
"    outline: 0;\n"
"	border: 1px solid gray;\n"
"\n"
"}\n"
"\n"
"\n"
"QListView::disabled \n"
"{\n"
"	background-color: #656565;\n"
"	color: #1b1b1b;\n"
"    border: 1px solid #656565;\n"
"\n"
"}\n"
"\n"
"\n"
"QListView::item \n"
"{\n"
""
                        "	background-color: #2d2d2d;\n"
"    padding: 1px;\n"
"\n"
"}\n"
"\n"
"\n"
"QListView::item:alternate \n"
"{\n"
"    background-color: #3a3a3a;\n"
"\n"
"}\n"
"\n"
"\n"
"QListView::item:selected \n"
"{\n"
"	background-color: #b78620;\n"
"	border: 1px solid #b78620;\n"
"	color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"QListView::item:selected:!active \n"
"{\n"
"	background-color: #b78620;\n"
"	border: 1px solid #b78620;\n"
"	color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"QListView::item:selected:active \n"
"{\n"
"	background-color: #b78620;\n"
"	border: 1px solid #b78620;\n"
"	color: #fff;\n"
"\n"
"}\n"
"\n"
"\n"
"QListView::item:hover {\n"
"    background-color: #262626;\n"
"    border: none;\n"
"    color: white;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QCheckBox-----*/\n"
"QCheckBox\n"
"{\n"
"	background-color: transparent;\n"
"    color: lightgray;\n"
"	border: none;\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator\n"
"{\n"
"    background-color: #323232;\n"
"    border: 1px solid darkgray;\n"
"    width: 12px;\n"
"    height: 12px;\n"
"\n"
""
                        "}\n"
"\n"
"\n"
"QCheckBox::indicator:checked\n"
"{\n"
"    image:url(\"./images/check.png\");\n"
"	background-color: #b78620;\n"
"    border: 1px solid #3a546e;\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:unchecked:hover\n"
"{\n"
"	border: 1px solid #b78620; \n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::disabled\n"
"{\n"
"	color: #656565;\n"
"\n"
"}\n"
"\n"
"\n"
"QCheckBox::indicator:disabled\n"
"{\n"
"	background-color: #656565;\n"
"	color: #656565;\n"
"    border: 1px solid #656565;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QRadioButton-----*/\n"
"QRadioButton \n"
"{\n"
"	color: lightgray;\n"
"	background-color: transparent;\n"
"\n"
"}\n"
"\n"
"\n"
"QRadioButton::indicator::unchecked:hover \n"
"{\n"
"	background-color: lightgray;\n"
"	border: 2px solid #b78620;\n"
"	border-radius: 6px;\n"
"}\n"
"\n"
"\n"
"QRadioButton::indicator::checked \n"
"{\n"
"	border: 2px solid #b78620;\n"
"	border-radius: 6px;\n"
"	background-color: rgba(183,134,32,20%);  \n"
"	width: 9px; \n"
"	height: 9px; \n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QSlider---"
                        "--*/\n"
"QSlider::groove:horizontal \n"
"{\n"
"	background-color: transparent;\n"
"	height: 3px;\n"
"\n"
"}\n"
"\n"
"\n"
"QSlider::sub-page:horizontal \n"
"{\n"
"	background-color: #b78620;\n"
"\n"
"}\n"
"\n"
"\n"
"QSlider::add-page:horizontal \n"
"{\n"
"	background-color: #131313;\n"
"\n"
"}\n"
"\n"
"\n"
"QSlider::handle:horizontal \n"
"{\n"
"	background-color: #b78620;\n"
"	width: 14px;\n"
"	margin-top: -6px;\n"
"	margin-bottom: -6px;\n"
"	border-radius: 6px;\n"
"\n"
"}\n"
"\n"
"\n"
"QSlider::handle:horizontal:hover \n"
"{\n"
"	background-color: #d89e25;\n"
"	border-radius: 6px;\n"
"\n"
"}\n"
"\n"
"\n"
"QSlider::sub-page:horizontal:disabled \n"
"{\n"
"	background-color: #bbb;\n"
"	border-color: #999;\n"
"\n"
"}\n"
"\n"
"\n"
"QSlider::add-page:horizontal:disabled \n"
"{\n"
"	background-color: #eee;\n"
"	border-color: #999;\n"
"\n"
"}\n"
"\n"
"\n"
"QSlider::handle:horizontal:disabled \n"
"{\n"
"	background-color: #eee;\n"
"	border: 1px solid #aaa;\n"
"	border-radius: 3px;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QScrol"
                        "lBar-----*/\n"
"QScrollBar:horizontal\n"
"{\n"
"    border: 1px solid #222222;\n"
"    background-color: #3d3d3d;\n"
"    height: 15px;\n"
"    margin: 0px 16px 0 16px;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::handle:horizontal\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(97, 97, 97, 255),stop:1 rgba(90, 90, 90, 255));\n"
"	border: 1px solid #2d2d2d;\n"
"    min-height: 20px;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::add-line:horizontal\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(97, 97, 97, 255),stop:1 rgba(90, 90, 90, 255));\n"
"	border: 1px solid #2d2d2d;\n"
"    width: 15px;\n"
"    subcontrol-position: right;\n"
"    subcontrol-origin: margin;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::sub-line:horizontal\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(97, 97, 97, 255),stop:1 rgba(90, 90, 90, 255));\n"
"	border: 1px solid #2d2d2d;\n"
"    width: 15px;\n"
"    subcontrol-p"
                        "osition: left;\n"
"    subcontrol-origin: margin;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::right-arrow:horizontal\n"
"{\n"
"    image: url(:/images/arrow-right.png);\n"
"    width: 6px;\n"
"    height: 6px;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::left-arrow:horizontal\n"
"{\n"
"    image: url(:/images/arrow-left.png);\n"
"    width: 6px;\n"
"    height: 6px;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal\n"
"{\n"
"    background: none;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar:vertical\n"
"{\n"
"    background-color: #3d3d3d;\n"
"    width: 16px;\n"
"	border: 1px solid #2d2d2d;\n"
"    margin: 16px 0px 16px 0px;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::handle:vertical\n"
"{\n"
"    background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(97, 97, 97, 255),stop:1 rgba(90, 90, 90, 255));\n"
"	border: 1px solid #2d2d2d;\n"
"    min-height: 20px;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::add-line:vertical\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1"
                        ":1, y1:0, x2:1, y2:1, stop:0 rgba(97, 97, 97, 255),stop:1 rgba(90, 90, 90, 255));\n"
"	border: 1px solid #2d2d2d;\n"
"    height: 15px;\n"
"    subcontrol-position: bottom;\n"
"    subcontrol-origin: margin;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::sub-line:vertical\n"
"{\n"
"	background-color: qlineargradient(spread:repeat, x1:1, y1:0, x2:1, y2:1, stop:0 rgba(97, 97, 97, 255),stop:1 rgba(90, 90, 90, 255));\n"
"	border: 1px solid #2d2d2d;\n"
"    height: 15px;\n"
"    subcontrol-position: top;\n"
"    subcontrol-origin: margin;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::up-arrow:vertical\n"
"{\n"
"    image: url(:/images/arrow-up.png);\n"
"    width: 6px;\n"
"    height: 6px;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::down-arrow:vertical\n"
"{\n"
"    image: url(:/images/arrow-down.png);\n"
"    width: 6px;\n"
"    height: 6px;\n"
"\n"
"}\n"
"\n"
"\n"
"QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical\n"
"{\n"
"    background: none;\n"
"\n"
"}\n"
"\n"
"\n"
"/*-----QProgressBar-----*/\n"
"QProgressBar\n"
"{\n"
"    "
                        "border: 1px solid #666666;\n"
"    text-align: center;\n"
"	color: #000;\n"
"	font-weight: bold;\n"
"\n"
"}\n"
"\n"
"\n"
"QProgressBar::chunk\n"
"{\n"
"    background-color: #b78620;\n"
"    width: 30px;\n"
"    margin: 0.5px;\n"
"\n"
"}"));
        actionExit = new QAction(Assets);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionExport_obj = new QAction(Assets);
        actionExport_obj->setObjectName(QString::fromUtf8("actionExport_obj"));
        actionExport_ma = new QAction(Assets);
        actionExport_ma->setObjectName(QString::fromUtf8("actionExport_ma"));
        actionLoad_ma = new QAction(Assets);
        actionLoad_ma->setObjectName(QString::fromUtf8("actionLoad_ma"));
        actionGL = new QAction(Assets);
        actionGL->setObjectName(QString::fromUtf8("actionGL"));
        actionVectorized = new QAction(Assets);
        actionVectorized->setObjectName(QString::fromUtf8("actionVectorized"));
        actionExport_svg = new QAction(Assets);
        actionExport_svg->setObjectName(QString::fromUtf8("actionExport_svg"));
        actionExport_mi = new QAction(Assets);
        actionExport_mi->setObjectName(QString::fromUtf8("actionExport_mi"));
        action128x128 = new QAction(Assets);
        action128x128->setObjectName(QString::fromUtf8("action128x128"));
        action256x256 = new QAction(Assets);
        action256x256->setObjectName(QString::fromUtf8("action256x256"));
        action512x512 = new QAction(Assets);
        action512x512->setObjectName(QString::fromUtf8("action512x512"));
        action1024x1024 = new QAction(Assets);
        action1024x1024->setObjectName(QString::fromUtf8("action1024x1024"));
        action2048x2048 = new QAction(Assets);
        action2048x2048->setObjectName(QString::fromUtf8("action2048x2048"));
        centralwidget = new QWidget(Assets);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        hboxLayout = new QHBoxLayout(centralwidget);
        hboxLayout->setSpacing(0);
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        Objects_widget = new QWidget(centralwidget);
        Objects_widget->setObjectName(QString::fromUtf8("Objects_widget"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(Objects_widget->sizePolicy().hasHeightForWidth());
        Objects_widget->setSizePolicy(sizePolicy);
        Objects_widget->setMinimumSize(QSize(320, 230));
        verticalLayout = new QVBoxLayout(Objects_widget);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        groupBox = new QGroupBox(Objects_widget);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setMinimumSize(QSize(0, 120));
        gridLayout_2 = new QGridLayout(groupBox);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        btn_addDetails = new QPushButton(groupBox);
        btn_addDetails->setObjectName(QString::fromUtf8("btn_addDetails"));

        gridLayout_2->addWidget(btn_addDetails, 2, 1, 1, 1);

        btn_loadShader = new QPushButton(groupBox);
        btn_loadShader->setObjectName(QString::fromUtf8("btn_loadShader"));

        gridLayout_2->addWidget(btn_loadShader, 2, 0, 1, 1);

        openGaussianFile = new QPushButton(groupBox);
        openGaussianFile->setObjectName(QString::fromUtf8("openGaussianFile"));

        gridLayout_2->addWidget(openGaussianFile, 0, 1, 1, 1);

        btn_saveGaussian = new QPushButton(groupBox);
        btn_saveGaussian->setObjectName(QString::fromUtf8("btn_saveGaussian"));

        gridLayout_2->addWidget(btn_saveGaussian, 1, 1, 1, 1);

        btn_exporthighres = new QPushButton(groupBox);
        btn_exporthighres->setObjectName(QString::fromUtf8("btn_exporthighres"));

        gridLayout_2->addWidget(btn_exporthighres, 4, 0, 1, 1);

        btn_openTerrain = new QPushButton(groupBox);
        btn_openTerrain->setObjectName(QString::fromUtf8("btn_openTerrain"));

        gridLayout_2->addWidget(btn_openTerrain, 1, 0, 1, 1);

        btn_clear = new QPushButton(groupBox);
        btn_clear->setObjectName(QString::fromUtf8("btn_clear"));

        gridLayout_2->addWidget(btn_clear, 0, 0, 1, 1);

        verticalSpacer = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Minimum);

        gridLayout_2->addItem(verticalSpacer, 3, 0, 1, 1);

        combo_exporthighres = new QComboBox(groupBox);
        combo_exporthighres->addItem(QString());
        combo_exporthighres->addItem(QString());
        combo_exporthighres->addItem(QString());
        combo_exporthighres->addItem(QString());
        combo_exporthighres->addItem(QString());
        combo_exporthighres->addItem(QString());
        combo_exporthighres->setObjectName(QString::fromUtf8("combo_exporthighres"));

        gridLayout_2->addWidget(combo_exporthighres, 4, 1, 1, 1);


        verticalLayout->addWidget(groupBox);

        groupBox_3 = new QGroupBox(Objects_widget);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        groupBox_3->setMinimumSize(QSize(0, 150));
        check_showInfluence = new QCheckBox(groupBox_3);
        check_showInfluence->setObjectName(QString::fromUtf8("check_showInfluence"));
        check_showInfluence->setGeometry(QRect(10, 70, 201, 23));
        slider_nbgaussians = new QSlider(groupBox_3);
        slider_nbgaussians->setObjectName(QString::fromUtf8("slider_nbgaussians"));
        slider_nbgaussians->setGeometry(QRect(10, 50, 271, 20));
        slider_nbgaussians->setMinimum(1);
        slider_nbgaussians->setMaximum(150);
        slider_nbgaussians->setValue(1);
        slider_nbgaussians->setOrientation(Qt::Horizontal);
        slider_nbgaussians->setTickPosition(QSlider::NoTicks);
        slider_nbgaussians->setTickInterval(0);
        label_nbgaussians = new QLabel(groupBox_3);
        label_nbgaussians->setObjectName(QString::fromUtf8("label_nbgaussians"));
        label_nbgaussians->setGeometry(QRect(200, 30, 81, 20));
        label_3 = new QLabel(groupBox_3);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(10, 30, 71, 16));
        label_9 = new QLabel(groupBox_3);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        label_9->setGeometry(QRect(10, 100, 91, 16));
        label_7 = new QLabel(groupBox_3);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(150, 100, 91, 16));
        slider_amplitude = new QSlider(groupBox_3);
        slider_amplitude->setObjectName(QString::fromUtf8("slider_amplitude"));
        slider_amplitude->setGeometry(QRect(10, 120, 131, 20));
        slider_amplitude->setMinimum(10);
        slider_amplitude->setMaximum(800);
        slider_amplitude->setValue(200);
        slider_amplitude->setOrientation(Qt::Horizontal);
        slider_amplitude->setTickPosition(QSlider::NoTicks);
        slider_amplitude->setTickInterval(0);
        slider_noiseLevel = new QSlider(groupBox_3);
        slider_noiseLevel->setObjectName(QString::fromUtf8("slider_noiseLevel"));
        slider_noiseLevel->setGeometry(QRect(150, 120, 131, 20));
        slider_noiseLevel->setMinimum(0);
        slider_noiseLevel->setMaximum(1000);
        slider_noiseLevel->setValue(500);
        slider_noiseLevel->setOrientation(Qt::Horizontal);
        slider_noiseLevel->setTickPosition(QSlider::NoTicks);
        slider_noiseLevel->setTickInterval(0);

        verticalLayout->addWidget(groupBox_3);

        groupBox_2 = new QGroupBox(Objects_widget);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        groupBox_2->setMinimumSize(QSize(0, 90));
        gridLayout_3 = new QGridLayout(groupBox_2);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        btn_graphTool = new QPushButton(groupBox_2);
        btn_graphTool->setObjectName(QString::fromUtf8("btn_graphTool"));

        gridLayout_3->addWidget(btn_graphTool, 1, 1, 1, 1);

        btn_handTool = new QPushButton(groupBox_2);
        btn_handTool->setObjectName(QString::fromUtf8("btn_handTool"));

        gridLayout_3->addWidget(btn_handTool, 0, 0, 1, 1);

        btn_eraseTool = new QPushButton(groupBox_2);
        btn_eraseTool->setObjectName(QString::fromUtf8("btn_eraseTool"));

        gridLayout_3->addWidget(btn_eraseTool, 0, 1, 1, 1);

        btn_moveTool = new QPushButton(groupBox_2);
        btn_moveTool->setObjectName(QString::fromUtf8("btn_moveTool"));

        gridLayout_3->addWidget(btn_moveTool, 1, 0, 1, 1);


        verticalLayout->addWidget(groupBox_2);

        toolOptionsWidget = new QStackedWidget(Objects_widget);
        toolOptionsWidget->setObjectName(QString::fromUtf8("toolOptionsWidget"));
        page_2 = new QWidget();
        page_2->setObjectName(QString::fromUtf8("page_2"));
        toolOptionsWidget->addWidget(page_2);
        page_3 = new QWidget();
        page_3->setObjectName(QString::fromUtf8("page_3"));
        groupBox_brush = new QGroupBox(page_3);
        groupBox_brush->setObjectName(QString::fromUtf8("groupBox_brush"));
        groupBox_brush->setEnabled(true);
        groupBox_brush->setGeometry(QRect(0, 10, 301, 75));
        groupBox_brush->setMinimumSize(QSize(0, 75));
        slider_brushthreshold = new QSlider(groupBox_brush);
        slider_brushthreshold->setObjectName(QString::fromUtf8("slider_brushthreshold"));
        slider_brushthreshold->setGeometry(QRect(110, 50, 181, 20));
        slider_brushthreshold->setMinimum(0);
        slider_brushthreshold->setMaximum(200);
        slider_brushthreshold->setValue(20);
        slider_brushthreshold->setOrientation(Qt::Horizontal);
        slider_brushthreshold->setTickPosition(QSlider::NoTicks);
        slider_brushthreshold->setTickInterval(0);
        label_4 = new QLabel(groupBox_brush);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(110, 30, 91, 16));
        btn_applyTool = new QPushButton(groupBox_brush);
        btn_applyTool->setObjectName(QString::fromUtf8("btn_applyTool"));
        btn_applyTool->setGeometry(QRect(10, 40, 92, 28));
        groupBox_4 = new QGroupBox(page_3);
        groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
        groupBox_4->setGeometry(QRect(0, 90, 301, 191));
        list_brushes = new QListWidget(groupBox_4);
        list_brushes->setObjectName(QString::fromUtf8("list_brushes"));
        list_brushes->setGeometry(QRect(150, 40, 141, 141));
        btn_openBrush = new QPushButton(groupBox_4);
        btn_openBrush->setObjectName(QString::fromUtf8("btn_openBrush"));
        btn_openBrush->setGeometry(QRect(10, 70, 131, 28));
        btn_saveBrush = new QPushButton(groupBox_4);
        btn_saveBrush->setObjectName(QString::fromUtf8("btn_saveBrush"));
        btn_saveBrush->setGeometry(QRect(10, 40, 131, 28));
        toolOptionsWidget->addWidget(page_3);
        page_4 = new QWidget();
        page_4->setObjectName(QString::fromUtf8("page_4"));
        groupBox_graph = new QGroupBox(page_4);
        groupBox_graph->setObjectName(QString::fromUtf8("groupBox_graph"));
        groupBox_graph->setEnabled(true);
        groupBox_graph->setGeometry(QRect(0, 10, 301, 111));
        slider_depthGraph = new QSlider(groupBox_graph);
        slider_depthGraph->setObjectName(QString::fromUtf8("slider_depthGraph"));
        slider_depthGraph->setGeometry(QRect(10, 80, 131, 16));
        slider_depthGraph->setMinimum(0);
        slider_depthGraph->setMaximum(60);
        slider_depthGraph->setValue(1);
        slider_depthGraph->setOrientation(Qt::Horizontal);
        slider_depthGraph->setTickPosition(QSlider::NoTicks);
        slider_depthGraph->setTickInterval(0);
        label_5 = new QLabel(groupBox_graph);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(10, 60, 91, 16));
        label_8 = new QLabel(groupBox_graph);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(150, 60, 91, 16));
        slider_blendThresholdGraph = new QSlider(groupBox_graph);
        slider_blendThresholdGraph->setObjectName(QString::fromUtf8("slider_blendThresholdGraph"));
        slider_blendThresholdGraph->setGeometry(QRect(150, 80, 131, 16));
        slider_blendThresholdGraph->setMinimum(0);
        slider_blendThresholdGraph->setMaximum(100);
        slider_blendThresholdGraph->setValue(100);
        slider_blendThresholdGraph->setOrientation(Qt::Horizontal);
        slider_blendThresholdGraph->setTickPosition(QSlider::NoTicks);
        slider_blendThresholdGraph->setTickInterval(0);
        check_influence = new QCheckBox(groupBox_graph);
        check_influence->setObjectName(QString::fromUtf8("check_influence"));
        check_influence->setGeometry(QRect(10, 30, 161, 19));
        check_influence->setChecked(false);
        toolOptionsWidget->addWidget(page_4);
        page = new QWidget();
        page->setObjectName(QString::fromUtf8("page"));
        groupBox_edit = new QGroupBox(page);
        groupBox_edit->setObjectName(QString::fromUtf8("groupBox_edit"));
        groupBox_edit->setEnabled(true);
        groupBox_edit->setGeometry(QRect(0, 10, 301, 131));
        radio_erase = new QRadioButton(groupBox_edit);
        radio_erase->setObjectName(QString::fromUtf8("radio_erase"));
        radio_erase->setGeometry(QRect(10, 30, 110, 24));
        radio_amplitude = new QRadioButton(groupBox_edit);
        radio_amplitude->setObjectName(QString::fromUtf8("radio_amplitude"));
        radio_amplitude->setGeometry(QRect(140, 30, 110, 24));
        radio_warp = new QRadioButton(groupBox_edit);
        radio_warp->setObjectName(QString::fromUtf8("radio_warp"));
        radio_warp->setGeometry(QRect(10, 90, 110, 24));
        radio_move = new QRadioButton(groupBox_edit);
        radio_move->setObjectName(QString::fromUtf8("radio_move"));
        radio_move->setGeometry(QRect(140, 90, 110, 24));
        radio_amplitude_lr = new QRadioButton(groupBox_edit);
        radio_amplitude_lr->setObjectName(QString::fromUtf8("radio_amplitude_lr"));
        radio_amplitude_lr->setGeometry(QRect(10, 60, 110, 24));
        radio_amplitude_hr = new QRadioButton(groupBox_edit);
        radio_amplitude_hr->setObjectName(QString::fromUtf8("radio_amplitude_hr"));
        radio_amplitude_hr->setGeometry(QRect(140, 60, 110, 24));
        toolOptionsWidget->addWidget(page);

        verticalLayout->addWidget(toolOptionsWidget);

        check_saveLogs = new QCheckBox(Objects_widget);
        check_saveLogs->setObjectName(QString::fromUtf8("check_saveLogs"));
        check_saveLogs->setEnabled(true);
        check_saveLogs->setChecked(true);

        verticalLayout->addWidget(check_saveLogs);

        btn_curveTool = new QPushButton(Objects_widget);
        btn_curveTool->setObjectName(QString::fromUtf8("btn_curveTool"));
        btn_curveTool->setEnabled(true);

        verticalLayout->addWidget(btn_curveTool);

        btn_AdrienGraph = new QPushButton(Objects_widget);
        btn_AdrienGraph->setObjectName(QString::fromUtf8("btn_AdrienGraph"));
        btn_AdrienGraph->setEnabled(true);

        verticalLayout->addWidget(btn_AdrienGraph);

        btn_test = new QPushButton(Objects_widget);
        btn_test->setObjectName(QString::fromUtf8("btn_test"));
        btn_test->setEnabled(true);

        verticalLayout->addWidget(btn_test);


        hboxLayout->addWidget(Objects_widget);

        raytracingwidget = new QOpenGLWidget(centralwidget);
        raytracingwidget->setObjectName(QString::fromUtf8("raytracingwidget"));

        hboxLayout->addWidget(raytracingwidget);

        Assets->setCentralWidget(centralwidget);
        menuBar = new QMenuBar(Assets);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1462, 22));
        menuOptions = new QMenu(menuBar);
        menuOptions->setObjectName(QString::fromUtf8("menuOptions"));
        menuRender_resolution = new QMenu(menuOptions);
        menuRender_resolution->setObjectName(QString::fromUtf8("menuRender_resolution"));
        Assets->setMenuBar(menuBar);

        menuBar->addAction(menuOptions->menuAction());
        menuOptions->addAction(menuRender_resolution->menuAction());
        menuRender_resolution->addAction(action128x128);
        menuRender_resolution->addAction(action256x256);
        menuRender_resolution->addAction(action512x512);
        menuRender_resolution->addAction(action1024x1024);
        menuRender_resolution->addAction(action2048x2048);

        retranslateUi(Assets);

        toolOptionsWidget->setCurrentIndex(2);


        QMetaObject::connectSlotsByName(Assets);
    } // setupUi

    void retranslateUi(QMainWindow *Assets)
    {
        Assets->setWindowTitle(QCoreApplication::translate("Assets", "Primitives Editor", nullptr));
        actionExit->setText(QCoreApplication::translate("Assets", "Exit", nullptr));
        actionExport_obj->setText(QCoreApplication::translate("Assets", "OBJ File (.obj)", nullptr));
        actionExport_ma->setText(QCoreApplication::translate("Assets", "Maya File (.ma)", nullptr));
        actionLoad_ma->setText(QCoreApplication::translate("Assets", "Load", nullptr));
        actionGL->setText(QCoreApplication::translate("Assets", "GL", nullptr));
        actionVectorized->setText(QCoreApplication::translate("Assets", "Vectorized", nullptr));
        actionExport_svg->setText(QCoreApplication::translate("Assets", "SVG File (.svg)", nullptr));
        actionExport_mi->setText(QCoreApplication::translate("Assets", "Maya Instance (.mi)", nullptr));
        action128x128->setText(QCoreApplication::translate("Assets", "128x128", nullptr));
        action256x256->setText(QCoreApplication::translate("Assets", "256x256", nullptr));
        action512x512->setText(QCoreApplication::translate("Assets", "512x512", nullptr));
        action1024x1024->setText(QCoreApplication::translate("Assets", "1024x1024", nullptr));
        action2048x2048->setText(QCoreApplication::translate("Assets", "2048x2048", nullptr));
        groupBox->setTitle(QCoreApplication::translate("Assets", "Input / Output", nullptr));
        btn_addDetails->setText(QCoreApplication::translate("Assets", "Add Primitives Details", nullptr));
        btn_loadShader->setText(QCoreApplication::translate("Assets", "Reload shader", nullptr));
        openGaussianFile->setText(QCoreApplication::translate("Assets", "Load Primitives", nullptr));
        btn_saveGaussian->setText(QCoreApplication::translate("Assets", "Save Primitives", nullptr));
        btn_exporthighres->setText(QCoreApplication::translate("Assets", "Export high-res", nullptr));
        btn_openTerrain->setText(QCoreApplication::translate("Assets", "Load Heightfield", nullptr));
        btn_clear->setText(QCoreApplication::translate("Assets", "Empty Terrain", nullptr));
        combo_exporthighres->setItemText(0, QCoreApplication::translate("Assets", "128x128", nullptr));
        combo_exporthighres->setItemText(1, QCoreApplication::translate("Assets", "256x256", nullptr));
        combo_exporthighres->setItemText(2, QCoreApplication::translate("Assets", "512x512", nullptr));
        combo_exporthighres->setItemText(3, QCoreApplication::translate("Assets", "1024x1024", nullptr));
        combo_exporthighres->setItemText(4, QCoreApplication::translate("Assets", "2048x2048", nullptr));
        combo_exporthighres->setItemText(5, QCoreApplication::translate("Assets", "4096x4096", nullptr));

        combo_exporthighres->setCurrentText(QCoreApplication::translate("Assets", "128x128", nullptr));
        groupBox_3->setTitle(QCoreApplication::translate("Assets", "View parameters", nullptr));
        check_showInfluence->setText(QCoreApplication::translate("Assets", "Show influence of the primitives", nullptr));
        label_nbgaussians->setText(QCoreApplication::translate("Assets", "Placeholder", nullptr));
        label_3->setText(QCoreApplication::translate("Assets", "Nb primitives", nullptr));
        label_9->setText(QCoreApplication::translate("Assets", "Amplitude", nullptr));
        label_7->setText(QCoreApplication::translate("Assets", "Details level", nullptr));
        groupBox_2->setTitle(QCoreApplication::translate("Assets", "Tools", nullptr));
        btn_graphTool->setText(QCoreApplication::translate("Assets", "Crest/River Graph", nullptr));
        btn_handTool->setText(QCoreApplication::translate("Assets", "Hand", nullptr));
#if QT_CONFIG(tooltip)
        btn_eraseTool->setToolTip(QCoreApplication::translate("Assets", "Ctrl + Wheel : change radius size", nullptr));
#endif // QT_CONFIG(tooltip)
        btn_eraseTool->setText(QCoreApplication::translate("Assets", "Simple disc", nullptr));
        btn_moveTool->setText(QCoreApplication::translate("Assets", "Region", nullptr));
        groupBox_brush->setTitle(QCoreApplication::translate("Assets", "Region Parameters", nullptr));
        label_4->setText(QCoreApplication::translate("Assets", "Mask threshold", nullptr));
        btn_applyTool->setText(QCoreApplication::translate("Assets", "Apply tool", nullptr));
        groupBox_4->setTitle(QCoreApplication::translate("Assets", "Templates", nullptr));
        btn_openBrush->setText(QCoreApplication::translate("Assets", "Open", nullptr));
        btn_saveBrush->setText(QCoreApplication::translate("Assets", "Save", nullptr));
        groupBox_graph->setTitle(QCoreApplication::translate("Assets", "Graph", nullptr));
        label_5->setText(QCoreApplication::translate("Assets", "Depth", nullptr));
        label_8->setText(QCoreApplication::translate("Assets", "Blend threshold", nullptr));
        check_influence->setText(QCoreApplication::translate("Assets", "Influence region control", nullptr));
        groupBox_edit->setTitle(QCoreApplication::translate("Assets", "Edit", nullptr));
        radio_erase->setText(QCoreApplication::translate("Assets", "Erase", nullptr));
        radio_amplitude->setText(QCoreApplication::translate("Assets", "Amplitude", nullptr));
        radio_warp->setText(QCoreApplication::translate("Assets", "Warp", nullptr));
        radio_move->setText(QCoreApplication::translate("Assets", "Move", nullptr));
        radio_amplitude_lr->setText(QCoreApplication::translate("Assets", "Amplitude LR", nullptr));
        radio_amplitude_hr->setText(QCoreApplication::translate("Assets", "Amplitude HR", nullptr));
        check_saveLogs->setText(QCoreApplication::translate("Assets", "Save logs", nullptr));
        btn_curveTool->setText(QCoreApplication::translate("Assets", "Curve", nullptr));
        btn_AdrienGraph->setText(QCoreApplication::translate("Assets", "Adrien Graph", nullptr));
        btn_test->setText(QCoreApplication::translate("Assets", "Test button", nullptr));
        menuOptions->setTitle(QCoreApplication::translate("Assets", "Options", nullptr));
        menuRender_resolution->setTitle(QCoreApplication::translate("Assets", "Render resolution", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Assets: public Ui_Assets {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UIC_MAIN_H
