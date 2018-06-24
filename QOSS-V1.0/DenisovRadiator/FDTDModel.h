#pragma once
#ifndef FDTDMODEL_H
#define FDTDMODEL_H

/*
这是一个新的FDTD显示窗口，含有VTK的三维显示
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

#include <QVTKWidget.h>

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>

#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderWindow.h>
#include <vtkCaptionActor2D.h>
#include <vtkTextProperty.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkGPUVolumeRayCastMapper.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkImageShiftScale.h>
#include <vtkImageCast.h>
#include <vtkImageData.h>
#include "vtkImageActor.h"
#include <vtkVolumeProperty.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkInformation.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>

#include "../DenisovRadiator/CodeJin/qcustomplot.h"		
#include "../util/Constant_Var.h"	
#include "../DenisovRadiator/QTFDTDThread.h"

using namespace std;
using namespace userInterface;

class FDTDModel : public QWidget
{
	Q_OBJECT
public:
	FDTDModel(QWidget *parent = 0);

	~FDTDModel();

private slots:
	void OpenModelFile(void);
	void OpenExcFile(void);
	void transparencyChanged(int v);
	void SetFDTDprocess(void);
	void RunFDTDprocess(void);
public slots:
	void RecieveFreq(double _freq);
	void recieveInt(int _in);
	void recieveFloat(float _in);

private:
	void CreateButtons(void);
	void ButtonsOff(void);
	void ButtonsOn(void);
	void update(void);
	void updateFDTDField(void);

private:
	//vtk 显示窗口
	vtkSmartPointer<vtkOrientationMarkerWidget> widget1;
	QVTKWidget widget;
	vtkSmartPointer<vtkRenderWindowInteractor> interactor;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkAxesActor> axes;

	vtkSmartPointer<vtkVolume> volume;
	vtkSmartPointer<vtkSmartVolumeMapper> volumeMapper;

	vtkSmartPointer<vtkVolume> Fieldvolume;
	vtkSmartPointer<vtkSmartVolumeMapper> FieldvolumeMapper;
	bool Field3D;

	vtkSmartPointer<vtkImageActor> fieldInActor;
	bool excisopen;
	bool modelisopen;
	bool fdtdcalculated;
	bool aperturecalculated;

	QScrollBar * transparencyScrollBar;

	int min_gray, max_gray;
	QLabel * transparencyLabel;
	QLineEdit * transparencyLineEdit;
	float transparency;

	//操作项
	QGroupBox *Buttons;
	QPushButton *OpenModel;
	QPushButton *OpenExc;
	QPushButton *SetFreList;
	QPushButton *SetFDTD;
	QPushButton *CalculateFDTD;
	QPushButton *SetAperture;
	QPushButton *CalAperture;

	//计算参数
	double frequency;

	vector<double> frequecyList;
	//数据
	int Nx, Ny, Nz;
	int Shiftx, Shifty;
	double dx, dy, dz;
	double cx, cy, cz;
	vector<double> Eps;
	vector<double> FDTDField;
	int OMPNumber;
	
	//显示：
	QGroupBox *ParaBox;
	QLabel *FreqLabel;	QLineEdit *FreqEdit;
	QLabel *dxLabel;	QLineEdit *dxEdit;
	QLabel *dyLabel;	QLineEdit *dyEdit;
	QLabel *dzLabel;	QLineEdit *dzEdit;
	QLabel *NxLabel;	QLineEdit *NxEdit;
	QLabel *NyLabel;	QLineEdit *NyEdit;
	QLabel *NzLabel;	QLineEdit *NzEdit;
	QLabel *OMPLabel;	QLineEdit *OMPEdit;

	//进度
	QLabel *Message;
	QProgressBar *progressBar;
	//FDTD计算线程
	QTFDTDThread *FDTD;

	std::string modelfile;
	std::string excfile;

};

#endif // FDTDMODEL_H