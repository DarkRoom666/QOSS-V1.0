#include "../DenisovRadiator/CodeJin/SourceModeGenerationD.h"//这个干掉 这个已经用了，最好
#include "../DenisovRadiator/CodeJin/Mathematical_Functions_Jin.h"
#include "../Calculation/Mathematical_Functions.h"
#include "../DenisovRadiator/CodeJin/CRadiator.h"
#include "../DenisovRadiator/HighOrderRadiator.h"
#include "../DenisovRadiator/showDenisov.h"
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <cmath>
using namespace calculation;
using namespace std;

showDenisov::showDenisov(QWidget *parent)
	: QWidget(parent)
{
	calculated = false;
	//ui.setupUi(this);
	resize(1400, 800);
	setWindowTitle("DenisovDesigner");
	//这是画图板 好像是李赟单独写的类
	paintField = new QPaintField;

	CreatePlot();

	CreateButtons();

	//ParaTable = new QTableView();
	//ParaTable -> setMaximumWidth(300);
	//ParaTable ->column
	
	CreateBasicParas();

	CreateDesignParas();

	CreateCutParas();

	Message = new QLabel(tr("Message"));
	Message -> setMaximumWidth(300);
	Percent = new QLabel(tr("Percentage"));
	Percent -> setMaximumWidth(300);

	Denisov1 = new DenisovRadiator;
	//FDTD的离散间隔
	dis = 0.1e-3;
	N_spa = 10;

	connect(btn1, SIGNAL(clicked()), this, SLOT(on_btn1()));
	//信号和槽函数的参数必须具有相同的顺序和相同的类型
	//connect(btn2, SIGNAL(clicked()), this, SLOT(ChangeValue()));
	connect(this, SIGNAL(SendValue(double)), this, SLOT(ChangeValue(double)));
	//连接Denisov设计类的进度
	connect(Denisov1, SIGNAL(SendValue(double)), this, SLOT(ChangeValue(double)));
	//连接Denisov设计类的说明
	connect(Denisov1, SIGNAL(SendText(QString)), this, SLOT(ChangeText(QString)));
	//连接按钮关
	connect(Denisov1, SIGNAL(SendBtnClose()), this, SLOT(BtnClose()));
	//连接按钮开
	connect(Denisov1, SIGNAL(SendBtnOpen()), this, SLOT(BtnOpen()));
	connect(btn2, SIGNAL(clicked()), this, SLOT(on_btn2()));
	connect(btn3, SIGNAL(clicked()), this, SLOT(on_btn3()));
	connect(btn4, SIGNAL(clicked()), this, SLOT(on_btn4()));
	connect(btn5, SIGNAL(clicked()), this, SLOT(on_btn5()));
	connect(btn6, SIGNAL(clicked()), this, SLOT(on_btn6()));
	//连接Denisov计算得到的功率曲线
	connect(Denisov1, SIGNAL(SendCoefficients(double, double, double, double, int)), this, SLOT(RecieveCoefficients(double, double ,double, double, int)));
	//连接Denisov计算得到的辐射器切向场
	connect(Denisov1, SIGNAL(SendFieldUpdate()), this, SLOT(UpdateTanE()));
	//连接Denisov计算得到的辐射器表面电流
	connect(Denisov1, SIGNAL(SendCurrentUpdate()), this, SLOT(UpdateSurfaceJ()));
	//connect(Denisov1, SIGNAL(SendTangentialEField(std::vector<std::vector<std::complex<double>>>, std::vector<std::vector<std::complex<double>>>, int, int)),
	//	this, SLOT(RecieveTangentialEField(std::vector<std::vector<std::complex<double>>>, std::vector<std::vector<std::complex<double>>>, int, int)));


	QVBoxLayout * vBoxLayout = new QVBoxLayout;
	vBoxLayout->addWidget(BasicParas);
	vBoxLayout->addWidget(DesignParas);
	vBoxLayout->addWidget(CCParas);
	vBoxLayout->addWidget(Buttons);
	vBoxLayout->addWidget(Message);
	vBoxLayout->addWidget(Percent);

	//Logo
	DrawLogo = new QGroupBox;
	QLabel * imgLlabel;
	imgLlabel = new QLabel;
	imgLlabel->setPixmap(QPixmap("Qt/images/Denisov.png"));
	imgLlabel->setScaledContents(true);
	imgLlabel->setMaximumSize(250,150);
	QGridLayout * layout1 = new QGridLayout;
	layout1->addWidget(imgLlabel);
	//defGroupBox = new QGroupBox;
	DrawLogo->setLayout(layout1);
	DrawLogo->setTitle(tr("Definition methods"));

	//Horizontal
	//Curves
	QVBoxLayout * DrawZoneVC = new QVBoxLayout;
	DrawZoneVC->addWidget(PlotCurve);
	DrawZoneVC->addWidget(PlotPower);
	//Logo AndFields
	QVBoxLayout * DrawZoneVF = new QVBoxLayout;
	DrawZoneVF->addWidget(DrawLogo);
	DrawZoneVF->addWidget(paintField);
	QHBoxLayout * DrawZoneHU = new QHBoxLayout;
	DrawZoneHU->addLayout(DrawZoneVF);
	DrawZoneHU->addLayout(DrawZoneVC);
	QHBoxLayout * DrawZoneHD = new QHBoxLayout;
	DrawZoneHD->addWidget(DrawCurrent);
	
	//DrawZoneHD->addWidget(paintCurrent);
	//Vertrical
	QVBoxLayout * DrawZoneV = new QVBoxLayout;
	DrawZoneV->addLayout(DrawZoneHU);
	DrawZoneV->addLayout(DrawZoneHD);


	QHBoxLayout *hBoxLayout = new QHBoxLayout(this);
	hBoxLayout->addLayout(DrawZoneV);
	hBoxLayout->addLayout(vBoxLayout);

	/*//原2上2下布局
	QHBoxLayout * DrawZoneHU = new QHBoxLayout;
	DrawZoneHU->addWidget(paintField);
	DrawZoneHU->addWidget(PlotCurve);
	QHBoxLayout * DrawZoneHD = new QHBoxLayout;
	DrawZoneHD->addWidget(PlotPower);
	DrawZoneHD->addWidget(DrawCurrent);
	//Vertrical
	QVBoxLayout * DrawZoneV = new QVBoxLayout;
	DrawZoneV->addLayout(DrawZoneHU);
	DrawZoneV->addLayout(DrawZoneHD);
	
	QHBoxLayout *hBoxLayout = new QHBoxLayout(this);
	hBoxLayout->addLayout(DrawZoneV);
	hBoxLayout->addLayout(vBoxLayout);
	*/

}

showDenisov::~showDenisov()
{

}

void showDenisov::ChangeValue(double _in) {
	QString ss;
	ss.setNum(_in);
	Percent->setText(ss);
}

void showDenisov::ChangeText(QString _in) {
	Message->setText(_in);
}
//这个函数进行初始化
void showDenisov::on_btn1()
{//按按钮1，读取基本参数，绘制输入模式场图，并计算基本参数

	ReadBasicParas();
	emit SendFreq(frequency);
	//绘制输入模式的场分布
	updatePaint();
	paintField->setWindowTitle("instantaneous E-Field");
	paintField->show();

	this->Nx = 250; this->Ny = 250;
	SourceModeGenerationD SourceMode(2,1,2,m,n,frequency,radius,0,0,Nx);
	SourceMode.FieldCalculation_Circular();
	vector<vector<double>> Absdata(Ny, vector<double>(Nx, 0));
	vector<vector<complex<double>>> FieldDatax(Ny, vector<complex<double>>(Nx, 0));
	vector<vector<complex<double>>> FieldDatay(Ny, vector<complex<double>>(Nx, 0));
	SourceMode.GetEX(FieldDatax);
	SourceMode.GetEY(FieldDatay);
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			Absdata[j][i] = sqrt(FieldDatax[j][i].real()*FieldDatax[j][i].real() + FieldDatay[j][i].real()*FieldDatay[j][i].real());
		}
	}

	paintField->setNum(Ny, Nx);
	paintField->setData(Absdata);
	updatePaint();

	//完成绘图，下面更新相应的参数 lcut delbeta1, delbeta2;
	double temp1, temp2, temp3;
	SourceMode.GetCircularWaveguideProperty(temp1, temp2, temp3, lcut);
	//lcut = lcut*1e3;
	//CutHeightEdit->setEnabled(true);
	CutHeightEdit->setText(QString::number(lcut*1e3, 10, 4));
	//CutHeightEdit->setEnabled(false);
	delbeta1 = (kzmnTE(m - 1, n,lambda,radius) - kzmnTE(m, n,lambda,radius))*1.03;
	//delbeta1 = ( kzmnTE(m-1, n, lambda, radius) - kzmnTE(m+1, n, lambda, radius) ) / 2;
	delbeta2 = -((kzmnTE(m,n,lambda,radius) - kzmnTE(m+3,n-1,lambda,radius))*0.5 + (kzmnTE(m,n,lambda,radius) - kzmnTE(m-3,n+1,lambda,radius))*0.5);
	delbeta1 = delbeta1 / 1e3;	delbeta2 = delbeta2 / 1e3;
	delBeta1Edit->setText(QString::number(delbeta1, 10, 4));
	delBeta2Edit->setText(QString::number(delbeta2, 10, 4));

	PlotCurve->clearGraphs();
	SourceMode.~SourceModeGenerationD();
	
}
//这个函数设定计算参数
void showDenisov::on_btn2()
{	
	ReadDesignParas();
	//Pass parameters into the Denisov Class
	Denisov1->SetupDenisov(this->frequency, this->radius, this->m, this->n);
	Denisov1->SetupDenisovParas(this->delbeta1, this->delbeta2, this->l1, this->l2, this->zc1, this->zc2, this->ls1, this->ls2, this->lc1, this->lc2, this->Hz, this->dz, this->Nz, this->mag1, this->mag2);
	Denisov1->SetTangentialEField(this->Nx, this->Ny);
	if (this->Nz % 6 == 0) { this->Sparse = 6; }
	else if (this->Nz % 4 == 0) { this->Sparse = 4; }

	this->Nphi = 240;	this->Nheight = this->Nz / this->Sparse;
	Denisov1->SetSurfaceCurrent(this->Nphi,this->Sparse);
	this->SJ.resize(this->Nphi);
	for (int p = 0; p < this->Nphi; p++) {
		this->SJ[p].resize(this->Nheight);
	}

	// 绘制扰动曲线
	std::vector<double> sig1, sig2;
	sig1.resize(Nz);	sig2.resize(Nz);
	Denisov1->GetDenisovTurbulence(sig1, sig2);

	QVector <double> Tx, T1y, T2y;
	Tx.resize(Nz);	T1y.resize(Nz);	T2y.resize(Nz);
	for (int i = 0; i < this->Nz; i++) {
		Tx[i] = (i + 0.5)*dz;
		T1y[i] = sig1[i];
		T2y[i] = sig2[i];
	}

	for (int i = 0; i < this->Nz; i++) {
		Tx[i] = Tx[i] * 1.0e3;
		T1y[i] = T1y[i] * 1.0e3;
		T2y[i] = T2y[i] * 1.0e3;
	}
	
	sig1.resize(0);	sig2.resize(0);
	//绘制扰动曲线
	PlotCurve->clearGraphs();
	PlotCurve->addGraph();
	PlotCurve->graph(0)->setPen(QPen(QColor(255, 110, 40)));//Red
	PlotCurve->graph(0)->setData(Tx, T1y);
	PlotCurve->addGraph();
	PlotCurve->graph(1)->setPen(QPen(QColor(40, 110, 255)));//Blue
	PlotCurve->graph(1)->setData(Tx, T2y);
	PlotCurve->xAxis->setRangeUpper(Hz*1e3);
	this->PlotCurve->replot();


	//绘制功率分布曲线-初值
	PowerTotal.resize(Nz+1);
	PowerMain.resize(Nz+1);
	PowerNeighbor.resize(Nz+1);
	PowerCorner.resize(Nz+1);
	ZAxis.resize(Nz+1);
	for (int i = 0; i < Nz+1; i++) {
		PowerTotal[i] = 100;
		PowerMain[i] = 100;
		PowerNeighbor[i] = 0;
		PowerCorner[i] = 0;
		ZAxis[i] = i*dz*1.0e3;
	}
	showDenisov::DrawPowerRatio();
	/*
	PlotPower->clearGraphs();
	PlotPower->xAxis->setRangeUpper(this->Hz * 1.0e3);
	PlotPower->addGraph();
	PlotPower->graph(0)->setPen(QPen(QColor(0, 0, 0)));
	PlotPower->graph(0)->setData(ZAxis, PowerTotal);
	PlotPower->graph(0)->setName("Total Modes Power");
	PlotPower->addGraph();
	PlotPower->graph(1)->setPen(QPen(QColor(255, 110, 40)));
	PlotPower->graph(1)->setData(ZAxis, PowerMain);
	PlotPower->graph(1)->setName("Central Modes Power");
	PlotPower->addGraph();
	PlotPower->graph(2)->setPen(QPen(QColor(40, 110, 255)));
	PlotPower->graph(2)->setData(ZAxis, PowerNeighbor);
	PlotPower->graph(2)->setName("Neighbor Modes Power");
	PlotPower->addGraph();
	PlotPower->graph(3)->setPen(QPen(QColor(110, 40, 255)));
	PlotPower->graph(3)->setData(ZAxis, PowerCorner);
	PlotPower->graph(3)->setName("Corner Modes Power");
	PlotPower->legend->setVisible(true);
	
	this->PlotPower->replot();
	*/
	

}
//这个函数启动DenisovCMT计算
void showDenisov::on_btn3()
{

	//Denisov1->run();

	Denisov1->start();
	//if (Denisov1->isFinished()) {
	//	QString message = "Denisov Simulation is done";
	//	emit SendText(message);
	//}

}
//这个函数绘制切口曲线
void showDenisov::on_btn4() {
	DrawCut();
	if (calculated) {
		OutputExc();
	}

}

void showDenisov::DrawCut() {
	bool ok = false;
	//这里是设定切口的位置
	showDenisov::zcut = showDenisov::CutPosEdit->text().toDouble(&ok);
	showDenisov::zcut = showDenisov::zcut*1e-3;
	showDenisov::TotalHEdit->setText(QString::number((zcut + lcut)*1e3, 10, 4));
	showDenisov::phic = (showDenisov::CutAngEdit->text().toInt(&ok));//取整！
	QString Write;	 Write.setNum(showDenisov::phic);
	showDenisov::CutAngEdit->setText(Write);
	if (showDenisov::zcut > 0 && showDenisov::zcut < showDenisov::Hz && showDenisov::phic>0 && showDenisov::phic < 360.0) {
		// 这里是画线的代码
		QVector<double> xv(361 + 10), yv(361 + 10);
		for (int i = 0; i < phic; i++) {
			xv[i] = zcut + lcut - phic / 360.0 * lcut + lcut*i / 360.0;
			xv[i] = xv[i] * 1e3;
			yv[i] = i;
		}
		for (int i = phic; i < phic + 10; i++) {
			xv[i] = zcut + lcut - (i - phic) / 10.0*lcut;
			xv[i] = xv[i] * 1e3;
			yv[i] = phic;
		}
		for (int i = phic + 10; i < 361 + 10; i++) {
			xv[i] = zcut + (i - 10 - phic) / 360.0 * lcut;
			xv[i] = xv[i] * 1e3;
			yv[i] = i - 10;
		}
		// create graph and assign data to it:
		DrawCurrent->addGraph();
		DrawCurrent->graph(0)->setData(xv, yv);
		DrawCurrent->replot();
	}
}
//这个函数存储计算结果
void showDenisov::on_btn5() {
	//First Define a DataStruct That Contains All the Parameter and Result 
	//Result include: Coefficient Curve Reults; Current Distribution Results; Cut Field Results;
	FILE* fileout;
	fileout = fopen("./DenisovParameters.dat", "wb");
	fwrite(&frequency,sizeof(double),1,fileout);
	fwrite(&radius, sizeof(double), 1, fileout);
	fwrite(&m, sizeof(int), 1, fileout);
	fwrite(&n, sizeof(int), 1, fileout);
	fwrite(&lcut, sizeof(double), 1, fileout);
	fwrite(&delbeta1, sizeof(double), 1, fileout);
	fwrite(&l1, sizeof(int), 1, fileout);
	fwrite(&mag1, sizeof(double), 1, fileout);
	fwrite(&ls1, sizeof(double), 1, fileout);
	fwrite(&zc1, sizeof(double),1, fileout);
	fwrite(&lc1, sizeof(double), 1, fileout);
	fwrite(&delbeta2, sizeof(double), 1, fileout);
	fwrite(&l2, sizeof(int), 1, fileout);
	fwrite(&mag2, sizeof(double), 1, fileout);
	fwrite(&ls2, sizeof(double), 1, fileout);
	fwrite(&zc2, sizeof(double), 1, fileout);
	fwrite(&lc2, sizeof(double), 1, fileout);
	fwrite(&Hz, sizeof(double), 1, fileout);
	fwrite(&dz, sizeof(double), 1, fileout);
	fwrite(&Nz, sizeof(int), 1, fileout);
	fwrite(&Sparse, sizeof(int), 1, fileout);
	fwrite(&phic, sizeof(int), 1, fileout);
	fwrite(&zcut, sizeof(double), 1, fileout);
	//Power Curves	
	for (int i = 0; i < Nz + 1; i++) {
		fwrite(&ZAxis[i], sizeof(double), 1, fileout);
	}
	for (int i = 0; i < Nz + 1; i++) {
		fwrite(&PowerTotal[i], sizeof(double), 1, fileout);
	}
	for (int i = 0; i < Nz + 1; i++) {
		fwrite(&PowerMain[i], sizeof(double), 1, fileout);
	}
	for (int i = 0; i < Nz + 1; i++) {
		fwrite(&PowerNeighbor[i], sizeof(double), 1, fileout);
	}
	for (int i = 0; i < Nz + 1; i++) {
		fwrite(&PowerCorner[i], sizeof(double), 1, fileout);
	}
	//Current Distribution
	for (int i = 0; i < Nphi; i++) {
		for (int j = 0; j < Nheight; j++) {
			fwrite(&SJ[i][j], sizeof(double), 1, fileout);
		}
	}
	fclose(fileout);

	OutputModel();
	OutputLattice();


}
//这个函数读取计算结果
void showDenisov::on_btn6() {
	FILE* filein;
	filein = fopen("./DenisovParameters.dat", "rb");
	fread(&frequency, sizeof(double), 1, filein);
	showDenisov::lambda = 2.998e8 / showDenisov::frequency;
	fread(&radius, sizeof(double), 1, filein);
	fread(&m, sizeof(int), 1, filein);
	fread(&n, sizeof(int), 1, filein);
	fread(&lcut, sizeof(double), 1, filein);
	fread(&delbeta1, sizeof(double), 1, filein);
	fread(&l1, sizeof(int), 1, filein);
	fread(&mag1, sizeof(double), 1, filein);
	fread(&ls1, sizeof(double), 1, filein);
	fread(&zc1, sizeof(double), 1, filein);
	fread(&lc1, sizeof(double), 1, filein);
	fread(&delbeta2, sizeof(double), 1, filein);
	fread(&l2, sizeof(int), 1, filein);
	fread(&mag2, sizeof(double), 1, filein);
	fread(&ls2, sizeof(double), 1, filein);
	fread(&zc2, sizeof(double), 1, filein);
	fread(&lc2, sizeof(double), 1, filein);
	fread(&Hz, sizeof(double), 1, filein);
	fread(&dz, sizeof(double), 1, filein);
	fread(&Nz, sizeof(int), 1, filein);
	fread(&Sparse, sizeof(int), 1, filein);
	fread(&phic, sizeof(int), 1, filein);
	fread(&zcut, sizeof(double), 1, filein);
	showDenisov::WriteAllParas();
	on_btn1();
	on_btn2();
	//Power Curves	
	for (int i = 0; i < Nz + 1; i++) {
		fread(&ZAxis[i], sizeof(double), 1, filein);
	}
	for (int i = 0; i < Nz + 1; i++) {
		fread(&PowerTotal[i], sizeof(double), 1, filein);
	}
	for (int i = 0; i < Nz + 1; i++) {
		fread(&PowerMain[i], sizeof(double), 1, filein);
	}
	for (int i = 0; i < Nz + 1; i++) {
		fread(&PowerNeighbor[i], sizeof(double), 1, filein);
	}
	for (int i = 0; i < Nz + 1; i++) {
		fread(&PowerCorner[i], sizeof(double), 1, filein);
	}
	//Current Distribution
	for (int i = 0; i < Nphi; i++) {
		for (int j = 0; j < Nheight; j++) {
			fread(&SJ[i][j], sizeof(double), 1, filein);
		}
	}
	fclose(filein);
	showDenisov::DrawPowerRatio();
	showDenisov::DrawSurfaceJ();
	calculated = false;
	on_btn4();
}

void showDenisov::WriteAllParas() {
	//FreqEdit->setText();
	FreqEdit->setText(QString::number(frequency/1.0e9, 10, 2));
	RadiusEdit->setText(QString::number(radius*1.0e3,10,2));
	MEdit->setText(QString::number(m, 10, 0));
	NEdit->setText(QString::number(n, 10, 0));
	CutHeightEdit->setText(QString::number(lcut*1.0e3, 10, 4));
	delBeta1Edit->setText(QString::number(delbeta1*1.0e-3, 10, 4));
	l1Edit->setText(QString::number(l1, 10, 0));
	mag1Edit->setText(QString::number(mag1*1.0e3, 10, 4));
	lc1Edit->setText(QString::number(lc1*1.0e3, 10, 4));
	ls1Edit->setText(QString::number(ls1*1.0e3, 10, 4));
	zc1Edit->setText(QString::number(zc1*1.0e3, 10, 4));
	delBeta2Edit->setText(QString::number(delbeta2*1.0e-3, 10, 4));
	l2Edit->setText(QString::number(l2, 10, 0));
	mag2Edit->setText(QString::number(mag2*1.0e3, 10, 4));
	lc2Edit->setText(QString::number(lc2*1.0e3, 10, 4));
	ls2Edit->setText(QString::number(ls2*1.0e3, 10, 4));
	zc2Edit->setText(QString::number(zc2*1.0e3, 10, 4));
	ComHeightEdit->setText(QString::number(Hz*1.0e3, 10, 4));
	ComNsEdit->setText(QString::number(Nz, 10, 0));
	CutPosEdit->setText(QString::number(zcut*1.0e3, 10, 2));
	CutAngEdit->setText(QString::number(phic, 10, 0));
	TotalHEdit->setText(QString::number((zcut + lcut)*1.0e3, 10, 4));
	
}

void showDenisov::RecieveCoefficients(double _CoeTotal, double _CoeMain, double _CoeNeighbor, double _CoeCorner,int _nn) {

	PowerTotal[_nn + 1] = _CoeTotal*100.0;
	PowerMain[_nn + 1] = _CoeMain*100.0;
	PowerNeighbor[_nn + 1] = _CoeNeighbor*100.0;
	PowerCorner[_nn + 1] = _CoeCorner*100.0;
	if (_nn % 10 == 0) {
		showDenisov::DrawPowerRatio();
	}
}

void showDenisov::DrawPowerRatio() {
	PlotPower->clearGraphs();
	PlotPower->addGraph();
	PlotPower->graph(0)->setPen(QPen(QColor(0, 0, 0)));
	PlotPower->graph(0)->setData(ZAxis, PowerTotal);
	PlotPower->graph(0)->setName("Total Modes Power");
	PlotPower->addGraph();
	PlotPower->graph(1)->setPen(QPen(QColor(255, 110, 40)));
	PlotPower->graph(1)->setData(ZAxis, PowerMain);
	PlotPower->graph(1)->setName("Central Modes Power");
	PlotPower->addGraph();
	PlotPower->graph(2)->setPen(QPen(QColor(40, 110, 255)));
	PlotPower->graph(2)->setData(ZAxis, PowerNeighbor);
	PlotPower->graph(2)->setName("Neighbor Modes Power");
	PlotPower->addGraph();
	PlotPower->graph(3)->setPen(QPen(QColor(110, 40, 255)));
	PlotPower->graph(3)->setData(ZAxis, PowerCorner);
	PlotPower->graph(3)->setName("Corner Modes Power");
	PlotPower->legend->setVisible(false);
	this->PlotPower->replot();
}

void showDenisov::UpdateTanE() {
	vector<vector<double>> Absdata(Nx, vector<double>(Ny, 0));
	vector<vector<complex<double>>> _Ex(Nx, vector<complex<double>>(Ny, 0));
	vector<vector<complex<double>>> _Ey(Nx, vector<complex<double>>(Ny, 0));
	Denisov1->GetTangentialEField(_Ex, _Ey, this->Nx, this->Ny);
	for (int j = 0; j < Nx; j++) {
		for (int i = 0; i < Ny; i++) {
			Absdata[j][i] = sqrt(_Ex[j][i].real()*_Ex[j][i].real() + _Ey[j][i].real()*_Ey[j][i].real());
		}
	}

	paintField->setNum(Nx, Ny);
	paintField->setData(Absdata);
	updatePaint();
}

void showDenisov::UpdateSurfaceJ() {

	Denisov1->GetSurfaceCurrent(SJ,Nphi,Nheight);

	showDenisov::DrawSurfaceJ();

	
}

void showDenisov::DrawSurfaceJ() {
	//
	QCPColorMap *colorMap = new QCPColorMap(this->DrawCurrent->xAxis, this->DrawCurrent->yAxis);
	int nx = this->Nheight;
	int ny = this->Nphi;
	colorMap->data()->setSize(nx, ny);
	//colorMap->data()->setRange(QCPRange(0, 100.0), QCPRange(0, 360.0)); // and span the coordinate range -4..4 in both key (x) and value (y) dimensions
	colorMap->data()->setRange(QCPRange(0, this->Hz*1e3), QCPRange(0, 360.0));
	double x, y, z;
	for (int xIndex = 0; xIndex<nx; xIndex++) {
		for (int yIndex = 0; yIndex<ny; yIndex++) {
			colorMap->data()->cellToCoord(xIndex, yIndex, &x, &y);
			z = SJ[yIndex][xIndex];
			colorMap->data()->setCell(xIndex, yIndex, z);
		}
	}

	// add a color scale:
	QCPColorScale *colorScale = new QCPColorScale(this->DrawCurrent);
	this->DrawCurrent->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
	colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
	colorMap->setColorScale(colorScale); // associate the color map with the color scale
	colorScale->axis()->setLabel("Induced Surface Current");

	// set the color gradient of the color map to one of the presets:
	colorMap->setGradient(QCPColorGradient::gpPolar);

	colorMap->rescaleDataRange();
	DrawCurrent->rescaleAxes();
	//DrawCurrent->xAxis->setVisible(false);
	//DrawCurrent->yAxis->setVisible(false);
	DrawCurrent->replot();
	//DrawCurrent->setInteraction(QCP::)
	//DrawCurrent->plottableDoubleClick;

}

void showDenisov::RecieveTangentialEField(std::vector<std::vector<std::complex<double>>> _Ex, std::vector<std::vector<std::complex<double>>> _Ey, int _Nx, int _Ny) {
	
	showDenisov::Nx = _Nx; showDenisov::Ny = _Ny;
	vector<vector<double>> E_Paint;
	E_Paint.resize(Nx);
	for (int i = 0; i < Nx; i++) {
		E_Paint[i].resize(Ny);
	}
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			E_Paint[i][j] = sqrt(abs(_Ex[i][j])*abs(_Ex[i][j]) + abs(_Ey[i][j])*abs(_Ey[i][j]));
		}
	}

	paintField->setNum(Nx, Ny);
	paintField->setData(E_Paint);
	updatePaint();

}

void showDenisov::BtnClose() {
	btn1->setEnabled(false);
	btn2->setEnabled(false);
	btn3->setEnabled(false);
	btn4->setEnabled(false);
	btn5->setEnabled(false);
	btn6->setEnabled(false);
	calculated = true;
}

void showDenisov::BtnOpen() {
	btn1->setEnabled(true);
	btn2->setEnabled(true);
	btn3->setEnabled(true);
	btn4->setEnabled(true);
	btn5->setEnabled(true);
	btn6->setEnabled(true);
}

void showDenisov::updatePaint()
{
	//每次更新时设置好新的数据  调用这个函数就行
	//paintField->setData(data);
	paintField->update();
	//paintField->setScaledContents(true);
	//paintField->setMaximumSize(150, 150);

}

void showDenisov::CreateButtons() {

	btn1 = new QPushButton(tr("SetBasic"));	btn1->setMaximumWidth(100);
	btn2 = new QPushButton(tr("SetDesign"));	btn2->setMaximumWidth(100);
	btn3 = new QPushButton(tr("RunCMT"));	btn3->setMaximumWidth(100);
	btn4 = new QPushButton(tr("SetCUT"));	btn4->setMaximumWidth(100);
	btn5 = new QPushButton(tr("Save"));		btn5->setMaximumWidth(100);
	btn6 = new QPushButton(tr("Load"));		btn6->setMaximumWidth(100);

	Buttons = new QGroupBox(tr("Buttons")); Buttons->setMaximumWidth(300);
	QVBoxLayout *ButtonV = new QVBoxLayout(Buttons);
	QHBoxLayout *ButtonH1 = new QHBoxLayout;
	QHBoxLayout *ButtonH2 = new QHBoxLayout;
	ButtonH1->addWidget(btn1);
	ButtonH1->addWidget(btn2);
	ButtonH1->addWidget(btn3);
	ButtonH2->addWidget(btn4);
	ButtonH2->addWidget(btn5);
	ButtonH2->addWidget(btn6);
	ButtonV->addLayout(ButtonH1);
	ButtonV->addLayout(ButtonH2);
}

void showDenisov::CreateBasicParas() {

	//Create BasicPara Aera
	FreqLabel = new QLabel(tr("Freq(GHz)"));		FreqLabel->setMaximumWidth(150);
	FreqEdit = new QLineEdit(tr("140"));			FreqEdit->setMaximumWidth(150);

	RadiusLabel = new QLabel(tr("Radius(mm)"));		RadiusLabel->setMaximumWidth(150);
	RadiusEdit = new QLineEdit(tr("16.8"));			RadiusEdit->setMaximumWidth(150);

	MLabel = new QLabel(tr("Mode m"));				MLabel->setMaximumWidth(150);
	MEdit = new QLineEdit(tr("22"));				MEdit->setMaximumWidth(150);

	NLabel = new QLabel(tr("Mode n"));				NLabel->setMaximumWidth(150);
	NEdit = new QLineEdit(tr("6"));					NEdit->setMaximumWidth(150);

	CutHeightLabel = new QLabel(tr("Cut Height, mm"));				CutHeightLabel->setMaximumWidth(150);
	CutHeightEdit = new QLineEdit(tr("Lc"));						CutHeightEdit->setMaximumWidth(150);	CutHeightEdit->setEnabled(false);

	BasicParas = new QGroupBox(tr("Basic Parameters"));		BasicParas->setMaximumWidth(300);
	QVBoxLayout *BasicV1 = new QVBoxLayout;
	BasicV1->addSpacing(10);
	BasicV1->addWidget(FreqLabel);
	BasicV1->addWidget(RadiusLabel);
	BasicV1->addWidget(MLabel);
	BasicV1->addWidget(NLabel);
	BasicV1->addWidget(CutHeightLabel);
	BasicV1->addSpacing(10);
	QVBoxLayout *BasicV2 = new QVBoxLayout;
	BasicV2->addSpacing(10);
	BasicV2->addWidget(FreqEdit);
	BasicV2->addWidget(RadiusEdit);
	BasicV2->addWidget(MEdit);
	BasicV2->addWidget(NEdit);
	BasicV2->addWidget(CutHeightEdit);
	BasicV2->addSpacing(10);
	QHBoxLayout *BasicH = new QHBoxLayout(BasicParas);
	BasicH->addLayout(BasicV1);
	BasicH->addLayout(BasicV2);

}

void showDenisov::CreateDesignParas() {
	
	//Labels and Edit
	delBeta1Label = new QLabel(tr("delbeta1,rad/mm"));		delBeta1Label->setMaximumWidth(150);
	delBeta1Edit = new QLineEdit(tr("50"));		delBeta1Edit->setMaximumWidth(150);
	delBeta2Label = new QLabel(tr("delbeta2,rad/mm"));		delBeta2Label->setMaximumWidth(150);
	delBeta2Edit = new QLineEdit(tr("20"));		delBeta2Edit->setMaximumWidth(150);
	l1Label = new QLabel(tr("l1"));					l1Label->setMaximumWidth(150);
	l1Edit = new QLineEdit(tr("1"));				l1Edit->setMaximumWidth(150);	l1Edit->setEnabled(false);
	l2Label = new QLabel(tr("l2"));					l2Label->setMaximumWidth(150);
	l2Edit = new QLineEdit(tr("3"));				l2Edit->setMaximumWidth(150);	l2Edit->setEnabled(false);
	mag1Label = new QLabel(tr("mag1 mm"));			mag1Label->setMaximumWidth(150);
	mag1Edit = new QLineEdit(tr("0.058"));			mag1Edit->setMaximumWidth(150);
	mag2Label = new QLabel(tr("mag2 mm"));			mag2Label->setMaximumWidth(150);
	mag2Edit = new QLineEdit(tr("0.047"));			mag2Edit->setMaximumWidth(150);
	zc1Label = new QLabel(tr("Central Position, turb1, mm"));	zc1Label->setMaximumWidth(150);
	zc1Edit = new QLineEdit(tr("25"));							zc1Edit->setMaximumWidth(150);
	zc2Label = new QLabel(tr("Central Position, turb2, mm"));	zc2Label->setMaximumWidth(150);
	zc2Edit = new QLineEdit(tr("35"));							zc2Edit->setMaximumWidth(150);
	lc1Label = new QLabel(tr("Central Length, turb1, mm"));	lc1Label->setMaximumWidth(150);
	lc1Edit = new QLineEdit(tr("30"));						lc1Edit->setMaximumWidth(150);
	lc2Label = new QLabel(tr("Central Length, turb2, mm"));	lc2Label->setMaximumWidth(150);
	lc2Edit = new QLineEdit(tr("40"));						lc2Edit->setMaximumWidth(150);
	ls1Label = new QLabel(tr("ChangeOver Length, turb1, mm"));	ls1Label->setMaximumWidth(150);
	ls1Edit = new QLineEdit(tr("10"));							ls1Edit->setMaximumWidth(150);
	ls2Label = new QLabel(tr("ChangeOver Length, turb2, mm"));	ls2Label->setMaximumWidth(150);
	ls2Edit = new QLineEdit(tr("10"));							ls2Edit->setMaximumWidth(150);
	//Compute Domain setup
	ComHeightLabel = new QLabel(tr("Compute Height, mm"));			CutHeightLabel->setMaximumWidth(150);
	ComHeightEdit = new QLineEdit(tr("150"));						CutHeightEdit->setMaximumWidth(150);
	ComNsLabel = new QLabel(tr("Compute Nsampling"));				ComNsLabel->setMaximumWidth(150);
	ComNsEdit = new QLineEdit(tr("3000"));							ComNsEdit->setMaximumWidth(150);

	DesignParas = new QGroupBox(tr("Turbulence Parameters"));
	DesignParas->setMaximumWidth(300);
	QVBoxLayout *LabelV = new QVBoxLayout;		QVBoxLayout *EditV = new QVBoxLayout;
	LabelV->addWidget(delBeta1Label);			EditV->addWidget(delBeta1Edit);
	LabelV->addWidget(l1Label);					EditV->addWidget(l1Edit);
	LabelV->addWidget(mag1Label);				EditV->addWidget(mag1Edit);
	LabelV->addWidget(zc1Label);				EditV->addWidget(zc1Edit);
	LabelV->addWidget(lc1Label);				EditV->addWidget(lc1Edit);
	LabelV->addWidget(ls1Label);				EditV->addWidget(ls1Edit);
	LabelV->addWidget(delBeta2Label);			EditV->addWidget(delBeta2Edit);
	LabelV->addWidget(l2Label);					EditV->addWidget(l2Edit);
	LabelV->addWidget(mag2Label);				EditV->addWidget(mag2Edit);
	LabelV->addWidget(zc2Label);				EditV->addWidget(zc2Edit);
	LabelV->addWidget(lc2Label);				EditV->addWidget(lc2Edit);
	LabelV->addWidget(ls2Label);				EditV->addWidget(ls2Edit);
	LabelV->addWidget(ComHeightLabel);			EditV->addWidget(ComHeightEdit);
	LabelV->addWidget(ComNsLabel);				EditV->addWidget(ComNsEdit);
	QHBoxLayout *GroupH = new QHBoxLayout(DesignParas);
	GroupH->addLayout(LabelV);
	GroupH->addLayout(EditV);
}

void showDenisov::CreateCutParas() {
	
	CutPosLabel = new QLabel(tr("Cut Height Position, mm"));		CutPosLabel->setMaximumWidth(150);
	CutPosEdit = new QLineEdit(tr("Zcut"));							CutPosEdit->setMaximumWidth(150);
	CutAngLabel = new QLabel(tr("Cut Angular Position, deg"));		CutAngLabel->setMaximumWidth(150);
	CutAngEdit = new QLineEdit(tr("Phic"));							CutAngEdit->setMaximumWidth(150);
	TotalHLabel	= new QLabel(tr("Total Height, mm"));				TotalHLabel->setMaximumWidth(150);
	TotalHEdit = new QLineEdit(tr("Zcut+Lc"));						TotalHEdit->setMaximumWidth(150);		TotalHEdit->setEnabled(false);
	CCParas = new QGroupBox("Cut Paras");		CCParas->setMaximumWidth(300);
	QVBoxLayout *LabelV = new QVBoxLayout;		QVBoxLayout *EditV = new QVBoxLayout;
	LabelV->addWidget(CutPosLabel);				EditV->addWidget(CutPosEdit);
	LabelV->addWidget(CutAngLabel);				EditV->addWidget(CutAngEdit);
	LabelV->addWidget(TotalHLabel);				EditV->addWidget(TotalHEdit);
	QHBoxLayout *GroupH = new QHBoxLayout(CCParas);
	GroupH->addLayout(LabelV);
	GroupH->addLayout(EditV);

}

void showDenisov::CreatePlot() {
	
	//Create Turbulence Curve Plot
	PlotCurve = new QCustomPlot;

//
	PlotCurve->xAxis->setLabel("height, z, mm");
	PlotCurve->yAxis->setLabel("turbulence mag, mm");
	// 设置背景色
	//PlotCurve->setBackground(QColor(50, 50, 50));
	PlotCurve->xAxis->setRange(0, 100);
	PlotCurve->yAxis->setRange(0, 0.1);
	PlotCurve->setMaximumHeight(200);
	PlotCurve->setMaximumWidth(400);

	//Create PowerDistribution Curve Plot
	PlotPower = new QCustomPlot;
	PlotPower->xAxis->setLabel("height,z,mm");
	PlotPower->yAxis->setLabel("Power Percentage, %");
	PlotPower->xAxis->setRange(0, 100);
	PlotPower->yAxis->setRange(0, 110);
	PlotPower->setMaximumHeight(200);
	PlotPower->setMaximumWidth(400);

	DrawCurrent = new QCustomPlot;
	DrawCurrent->xAxis->setLabel("height,z,mm");
	DrawCurrent->yAxis->setLabel("Angular,phi,deg");
	DrawCurrent->xAxis->setRange(0, 100);
	DrawCurrent->yAxis->setRange(0, 360);
	DrawCurrent->setMaximumHeight(270);
	DrawCurrent->setMaximumWidth(900);

}

void showDenisov::ReadBasicParas() {
	bool ok = false;
	//amplitude = valStr.toDouble(&ok);
	//读取基本参数
	showDenisov::frequency = showDenisov::FreqEdit->text().toDouble(&ok);
	showDenisov::frequency = showDenisov::frequency*1.0e9;	
	showDenisov::lambda = C_Speed / showDenisov::frequency;
	showDenisov::radius = showDenisov::RadiusEdit->text().toDouble(&ok);
	showDenisov::radius = showDenisov::radius*1.0e-3;
	showDenisov::m = showDenisov::MEdit->text().toInt(&ok);
	showDenisov::n = showDenisov::NEdit->text().toInt(&ok);
}

void showDenisov::ReadDesignParas() {
	//Read Design Parameters
	bool ok = false;

	showDenisov::l1 = showDenisov::l1Edit->text().toInt(&ok);
	showDenisov::l2 = showDenisov::l2Edit->text().toInt(&ok);
	showDenisov::delbeta1 = showDenisov::delBeta1Edit->text().toDouble(&ok);	showDenisov::delbeta1 = showDenisov::delbeta1 / 1.0e-3;
	showDenisov::delbeta2 = showDenisov::delBeta2Edit->text().toDouble(&ok);	showDenisov::delbeta2 = showDenisov::delbeta2 / 1.0e-3;
	showDenisov::zc1 = showDenisov::zc1Edit->text().toDouble(&ok);	showDenisov::zc1 = showDenisov::zc1*1e-3;
	showDenisov::zc2 = showDenisov::zc2Edit->text().toDouble(&ok);	showDenisov::zc2 = showDenisov::zc2*1e-3;
	showDenisov::lc1 = showDenisov::lc1Edit->text().toDouble(&ok);	showDenisov::lc1 = showDenisov::lc1*1e-3;
	showDenisov::lc2 = showDenisov::lc2Edit->text().toDouble(&ok);	showDenisov::lc2 = showDenisov::lc2*1e-3;
	showDenisov::ls1 = showDenisov::ls1Edit->text().toDouble(&ok);	showDenisov::ls1 = showDenisov::ls1*1e-3;
	showDenisov::ls2 = showDenisov::ls2Edit->text().toDouble(&ok);	showDenisov::ls2 = showDenisov::ls2*1e-3;
	showDenisov::mag1 = showDenisov::mag1Edit->text().toDouble(&ok);	showDenisov::mag1 = showDenisov::mag1*1e-3;
	showDenisov::mag2 = showDenisov::mag2Edit->text().toDouble(&ok);	showDenisov::mag2 = showDenisov::mag2*1e-3;
	showDenisov::Hz = showDenisov::ComHeightEdit->text().toDouble(&ok);	showDenisov::Hz = showDenisov::Hz*1e-3;
	showDenisov::Nz = showDenisov::ComNsEdit->text().toInt(&ok);		//dz = dz*1e-3;
	showDenisov::dz = showDenisov::Hz/ showDenisov::Nz;
}

void showDenisov::OutputModel() {

	Nx = (radius * 2) / dis;
	Ny = (radius * 2) / dis;
	//vector<float> Eps;
	
	//Eps.resize(Nx*Ny*Nz);
	fstream modelout;
	modelout.open("./Cutlog.txt",ios::out);

	CRadiator DenisovCut(2,1,2,showDenisov::m,showDenisov::n,showDenisov::Nx,showDenisov::Ny,float(showDenisov::lambda),
		                 showDenisov::radius,showDenisov::radius,showDenisov::lcut,float(N_spa*dis),0,0,0,0,0);
	DenisovCut.SetLogfile(&modelout);
	DenisovCut.ResetRadiator_RoundH(2,float(showDenisov::lambda),float(showDenisov::radius),showDenisov::Nx,showDenisov::Ny,N_spa,float(showDenisov::lcut),float(N_spa*dis));
	DenisovCut.GenerateCellArray_RoundH(float(showDenisov::phic*Pi/180.0));

	modelout.close();
	
	emit SendFreq(frequency);
	
}

void showDenisov::OutputExc() {
	Nx = (radius * 2) / dis;
	Ny = (radius * 2) / dis;
	int index = (zcut - N_spa*dis) / dz;
	vector<vector<complex<double>>> Ex;
	vector<vector<complex<double>>> Ey;
	vector<vector<complex<double>>> Hx;
	vector<vector<complex<double>>> Hy;

	Ex.resize(Nx); for (int i = 0; i < Nx; i++) { Ex[i].resize(Ny); }
	Ey.resize(Nx); for (int i = 0; i < Nx; i++) { Ey[i].resize(Ny); }
	Hx.resize(Nx); for (int i = 0; i < Nx; i++) { Hx[i].resize(Ny); }
	Hy.resize(Nx); for (int i = 0; i < Nx; i++) { Hy[i].resize(Ny); }

	Denisov1->GetExcitationField(Ex,Ey,Hx,Hy,index,Nx);
	//写文件

	fstream Fieldout;
	Fieldout.open("./DenisovExc.txt", ios::out);
	Fieldout << 0 << " " << 0 << " " << 0 << " "//Trans
			 << 0 << " " << 0 << " " << 1 << " " << 0 << " "//Roatate
		     << Nx << " " << Ny << " " << dis <<" "<<dis<< endl;
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			Fieldout
				<< abs(Ex[i][j]) << " " << arg(Ex[i][j]) * 180 / Pi << " "
				<< abs(Ey[i][j]) << " " << arg(Ey[i][j]) * 180 / Pi << " "
				<< 0 << " " << 0 << " "
				<< abs(Hx[i][j]) << " " << arg(Hx[i][j]) * 180 / Pi << " "
				<< abs(Hy[i][j]) << " " << arg(Hy[i][j]) * 180 / Pi << " "
				<< 0 << " " << 0 << " " << endl;
		}
	}
	Fieldout.close();
}

void showDenisov::OutputLattice() {
	HighOrderRadiator HOR;
	HOR.SetRadiatorBasicParas(radius,zcut,phic);
	HOR.SetRadiatorTurbulenceParas(delbeta1, delbeta2, zc1, zc2, ls1, ls2, lc1, lc2, mag1, mag2);
	SourceModeGenerationD SourceMode(2, 1, 2, m, n, frequency, radius, 0, 0, 100);
	double rc;
	double phi;
	double theta;
	double templcut;
	SourceMode.GetCircularWaveguideProperty(phi,theta,rc,templcut);
	HOR.SetFirstMirrorParas(0.04,0.16,rc,phi,lcut);
	SourceMode.~SourceModeGenerationD();


	FILE *Latticeout;
	
	vector<vector<Vector3>> Face;
	HOR.GetFirstMirrorLattice(Face, 2.998e8 / frequency / 2);
	int N1 = Face.size();
	int N2 = Face[0].size();
	Latticeout = fopen("./FaceLattice.dat","wb");
	fwrite(&N1, sizeof(int), 1, Latticeout);
	fwrite(&N2, sizeof(int), 1, Latticeout);
	for (int i = 0; i < N1; i++) {
		for (int j = 0; j < N2; j++) {
			float x = float(Face[i][j].x);
			float y = float(Face[i][j].y);
			float z = float(Face[i][j].z);
			fwrite(&x, sizeof(float), 1, Latticeout);
			fwrite(&y, sizeof(float), 1, Latticeout);
			fwrite(&z, sizeof(float), 1, Latticeout);
		}
	}
	fclose(Latticeout);

	vector<Vector3> Radiator1;
	HOR.GetRadiatorLattice1(Radiator1, 2.998e8 / frequency / 2);
	N1 = Radiator1.size();
	Latticeout = fopen("./RadiatorLattice1.dat", "wb");
	fwrite(&N1, sizeof(int), 1, Latticeout);
	for (int i = 0; i < N1; i++) {
			float x = float(Radiator1[i].x);
			float y = float(Radiator1[i].y);
			float z = float(Radiator1[i].z);
			fwrite(&x, sizeof(float), 1, Latticeout);
			fwrite(&y, sizeof(float), 1, Latticeout);
			fwrite(&z, sizeof(float), 1, Latticeout);
	}
	fclose(Latticeout);

	vector<Vector3> Radiator2;
	HOR.GetRadiatorLattice2(Radiator2, 2.998e8 / frequency / 2);
	N1 = Radiator2.size();
	Latticeout = fopen("./RadiatorLattice2.dat", "wb");
	fwrite(&N1, sizeof(int), 1, Latticeout);
	for (int i = 0; i < N1; i++) {
		float x = float(Radiator2[i].x);
		float y = float(Radiator2[i].y);
		float z = float(Radiator2[i].z);
		fwrite(&x, sizeof(float), 1, Latticeout);
		fwrite(&y, sizeof(float), 1, Latticeout);
		fwrite(&z, sizeof(float), 1, Latticeout);
	}
	fclose(Latticeout);
	
}

#include "moc_showDenisov.cpp"
