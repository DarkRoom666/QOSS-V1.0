#ifndef SHOWDENISOV_H
#define SHOWDENISOV_H

/*
这是显示主界面的类-可根据需要进行编辑
注意：这个是一个显示的窗口 只保留显示逻辑，
注意：所有运算部分调用DLL！ DLL剥离
*/

//这个是公共的
#include <QtWidgets/QWidget>
#include <QGroupBox>
#include <QPushButton>
#include <QTableView>
#include <QLineEdit>
#include <QString>
#include <iostream> 
#include <sstream> 
#include <string>
#include <windows.h>
#include <Vector>
#include <complex>

//#include "ui_showDenisov.h"		//这个是QT生成的


#include "CodeJin/qcustomplot.h"		//这个单独放到一个文件夹里
#include "CodeJin/QPaintField.h"		//这个也单独放到一个文件夹里

#include "../DenisovRadiator/DenisovRadiator.h"	//这个要封装成一个DLL 把SourceModeGenerationD集成进去最好-

#include "../util/Constant_Var.h"	//这个干掉

using namespace std;

class showDenisov : public QWidget
{
	Q_OBJECT

public:
	showDenisov(QWidget *parent = 0);
	~showDenisov();

private slots :

	void on_btn1();	//
	void on_btn2();
	void on_btn3();
	void on_btn4();
	void on_btn5();//Save
	void on_btn6();//Load
	void ChangeValue(double _in);
	void ChangeText(QString _in);
	void BtnClose();
	void BtnOpen();
	void RecieveCoefficients(double _CoeTotal, double _CoeMain, double _CoeNeighbor, double _CoeCorner, int _nn);
	void UpdateTanE();
	void UpdateSurfaceJ();
	void RecieveTangentialEField(std::vector<std::vector<std::complex<double>>> _Ex,std::vector<std::vector<std::complex<double>>> _Ey, int _Nx, int _Ny);
	//void DrawCurrents();


private: signals:
	void SendValue(double _out);
	void SendText(QString _Message);

public: signals:
	void SendFreq(double _freq);

private:
	void updatePaint();
	void CreatePlot();
	void CreateButtons();
	
	void CreateBasicParas();
	void CreateDesignParas();
	void CreateCutParas();
	void ReadBasicParas();
	void ReadDesignParas();
	void ReadCutParas();
	void WriteAllParas();
	void DrawSurfaceJ();
	void DrawPowerRatio();
	void DrawCut();
	void OutputModel();
	void OutputExc();
	void OutputLattice();
	void OutputOutField();

private:
	//Ui::showDenisovClass ui;
	bool calculated;
	QPaintField *paintField;	//Plot Input Mode Field
	QPaintField *paintCurrent;	//Plot Surface Current Distribution

	QCustomPlot *PlotCurve;		//Plot Turbulence Curves
	QCustomPlot *PlotPower;		//Plot PowerDistribution Curves;
	QCustomPlot *DrawCurrent;	//Plot Induced Surface Currents;
	QGroupBox *DrawLogo;

	QLabel *Percent;
	QLabel *Message;
	//Button Interface
	QGroupBox *Buttons;
	QPushButton *btn1;	//SetBasic
	QPushButton *btn2;	//SetDesign
	QPushButton *btn3;	//RunCMT
	QPushButton *btn4;	//SetCut
	QPushButton *btn5;	//Save
	QPushButton *btn6;	//Load

	//Basic Parameters InterFace
	QGroupBox *BasicParas;
	QLabel *FreqLabel;		QLineEdit *FreqEdit;
	QLabel *RadiusLabel;	QLineEdit *RadiusEdit;
	QLabel *MLabel;			QLineEdit *MEdit;
	QLabel *NLabel;			QLineEdit *NEdit;
	QLabel *CutHeightLabel;	QLineEdit *CutHeightEdit;
	
	//Design Parameters InterFace
	QGroupBox *DesignParas;
	QLabel *delBeta1Label;	QLineEdit *delBeta1Edit;
	QLabel *delBeta2Label;	QLineEdit *delBeta2Edit;
	QLabel *l1Label;		QLineEdit *l1Edit;
	QLabel *l2Label;		QLineEdit *l2Edit;
	QLabel *mag1Label;		QLineEdit *mag1Edit;
	QLabel *mag2Label;		QLineEdit *mag2Edit;
	QLabel *lc1Label;		QLineEdit *lc1Edit;
	QLabel *lc2Label;		QLineEdit *lc2Edit;
	QLabel *ls1Label;		QLineEdit *ls1Edit;
	QLabel *ls2Label;		QLineEdit *ls2Edit;
	QLabel *zc1Label;		QLineEdit *zc1Edit;
	QLabel *zc2Label;		QLineEdit *zc2Edit;
	QLabel *ComHeightLabel;	QLineEdit *ComHeightEdit;
	QLabel *ComNsLabel;		QLineEdit *ComNsEdit;

	//Cut Parameters InterFace
	QGroupBox *CCParas;
	QLabel *CutPosLabel;	QLineEdit *CutPosEdit;
	QLabel *CutAngLabel;	QLineEdit *CutAngEdit;
	QLabel *TotalHLabel;	QLineEdit *TotalHEdit;

	QTableView *ParaTable;
	DenisovRadiator *Denisov1;//CMT threads


//basic parameters
	double frequency;	double lambda;
	double radius;
	//mode
	int m;
	int n;
	//for Drawing Mode Field
	int Nx, Ny;
	int N_spa;
	//for Drawing Surface Current
	int Nphi, Nheight, Sparse;
	//Cut Height
	double lcut;
	//Cut Position
	int phic;
	double zcut;
//design parameters
	double delbeta1;
	double delbeta2;
	int l1;
	int l2;
	double mag1;
	double mag2;
	double ls1;
	double ls2;
	double zc1;
	double zc2;
	double lc1;
	double lc2;
	double Hz;
	double dz;
	double dis;
	int Nz;
//Compute parasmeters

//Display Data;

//Power
	QVector<double> PowerTotal;
	QVector<double> PowerMain;
	QVector<double> PowerNeighbor;
	QVector<double> PowerCorner;
	QVector<double> ZAxis;

	vector<vector<double>> SJ;

	vector<double> EpsModel;

	
};

#endif // SHOWFDTD_H
