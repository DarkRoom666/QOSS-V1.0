#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include <QVTKWidget.h>
#include <QAction>
#include <QTreeWidget>
#include <QDockWidget>
#include <QLabel>
#include <QPushButton>
#include <QMenu>
#include <QMenuBar>
#include <QToolBar>
#include <QstatusBar>
#include <QTreeWidgetItem>
#include <QButtonGroup>
#include <QRadioButton>

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>

#include "MyData.h"

#include "Qt/include/ParaboloidWidget.h"
#include "Qt/include/RestrictionWidget.h"
#include "Qt/include/GaussianWidget.h"
#include "Qt/include/ApertureFieldWidget.h"
#include "Qt/include/CalculationWidget.h"
#include "Qt/include/ParabolicCylinderWidget.h"
#include "Qt/include/PlaneMirrorWidget.h"
#include "Qt/include/STLMirrorWidget.h"
#include "Qt/include/GraphTransWidget.h"
#include "Qt/include/EllipsoidWidget.h"
#include "Qt/include/PhsCorProgressDialog.h"
#include "Qt/include/FDTDProgressDialog.h"
#include "Qt/include/PVVAProgressDialog.h"

#include <Calculation/FDTDRadiator.h>

using namespace userInterface;
class mainWindow : public QMainWindow
{
	Q_OBJECT
	enum QVariantType
	{	
		FIELD = 10,		
	};

public:
	mainWindow(QWidget *parent = 0);
	~mainWindow();

private:
	void init();

	//------------------- ��������------------------- 
	void createActions();
	void createMenus(); // �˵�
	void createToolBars();
	void createStatusBar();

	void createTreeWidgetItem(); // ����tree
	void createRightMenu(); // �Ҽ��˵�
	void createDetails(); //zuojian

	void createProject();
	// 
	void isNeedSave();

	void updateVtk();

	void updateLight();

	// �ı�3D��ʾ����
	void showDetails(int);


private slots:

	// �˵���Ӧ����
	void openFile();
	void newFile();
	void viewInitFile();
	void setView(int);

	void on_isShowBox();
	void on_isTransparentMirror();
	void on_isShowMirror();

	void on_modifyingRestriction();
	void on_delRestriction();

	void on_modifyingMirror();
	void on_modifyParameters();

	void on_restriction();
	void toReceiveRestriction(int);

	void toApertureField(int);
	void createApertureField();

	void createGaussian();
	void toReceiveGaussian(int);

	void on_createParaboloid();
	void on_createParabolicCylinder();
	void on_createEllipsoid();
	void on_createPlaneMirror();
	void on_createSTLMirror();

	void on_saveSTL();

	//void toReceiveParabolicCylinder(int);
	void toReceiveMirror(int);
	void toReceiveMirrorType(int);

	void on_PVVA();
	void on_FDTD();
	void toReceiveFDTD();
	void loadFDTDField(); //����FDTD����õĳ�
	void toReceiveFDTDStop();

	void on_PhaseCor();
	void toReceivePhaseCor(Mirror*);
	void toReceivePVVAField(Field*);

	// ------------------- �Ҽ����� ----------------------------------
	void on_treeWidget_ContextMenuRequested(QPoint pos);// �Ҽ��˵�

	// ------------------- ������� ----------------------------------
	void on_treeWidget_leftPressed(QTreeWidgetItem *item, int column);
	void on_Details_FieldClicked(); //details Field 


private:
	vtkSmartPointer<vtkOrientationMarkerWidget> widget1;
	QVTKWidget widget; // vtk ��ʾ����
	vtkSmartPointer<vtkRenderWindowInteractor> interactor; // ����
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkAxesActor> axes; // ����

	QDockWidget * leftWidget; //����treeWidget�Ĵ���
	QTreeWidget * treeWidget;

	QDockWidget * detailsDockWidget;
	QWidget * detailsWidget;

	QLabel * locationLabel;  // ״̬����ǩ
	QPushButton * unitBtn;   // ��λ��ť

	//----------- Menu ----------------- 
	QMenu * fileMenu;  //�˵���
	QMenu * viewMenu;  //��ʾ��
	QMenu * eidtMenu;  //�༭��
	QMenu * ModelMenu;  //ģ����
	QMenu * SourceMenu;  //Դ��
	QMenu * CalMenu;  //������

	//�Ҽ��˵�
	QMenu *R_Tree_MirrorTypeMenu;
	QMenu *R_Tree_MirrorParMenu;
	QMenu *R_Tree_RestrictionMenu;
	QMenu *R_BlankMenu;

	//----------- ToolBar ------------------- 
	QToolBar * fileTool;   //������
	QToolBar * constructTool;  //ģ����

	//----------- Action ----------------- 
	//�ļ��˵���
	QAction * saveFileAction;
	QAction * openFileAction;
	QAction * newFileAction;
	//�ļ��˵�--view
	QAction * viewAction;  // �ӽ�
	QLabel * viewLabel;
	QComboBox * viewComboBox;

	QAction * isShowBoxAction;

	QAction * GaussianAction;     // ��˹��Դ
	QAction * ApertureFieldAction;     // ���ⳡԴ
	QAction * PVVAAction;     // ����pVVA
	QAction * FDTDAction;     // ����FDTDAction
	QAction * PhaseAction;    // ��λ����

	// �Ҽ�
	QAction * modifyingMirrorAction;
	QAction * restrictionAction;
	QAction * modifyingRestrictionAction;
	QAction * delRestrictionAction;

	QAction * modifyParametersAction;
	QAction * isShowMirrorAction;
	QAction * isTransparentAction;
	QAction * saveSTLction;     // ����ģ��

	QAction * loadFDTDAction; // ����FDTD�����Դ

	//------------- TreeWidgetItem------------------- 
	QTreeWidgetItem * definitionsTreeItem;
	QTreeWidgetItem * variablesTreeItem;
	QTreeWidgetItem * modelTreeItem;
	QTreeWidgetItem * geometryTreeItem;
	QTreeWidgetItem * sourceTreeItem;
	QTreeWidgetItem * soucreFieldTreeItem;
	QTreeWidgetItem * lightTreeItem;
	QTreeWidgetItem * fieldTreeItem;

	QTreeWidgetItem * rightSelectItem;
	QTreeWidgetItem * leftSelectItem; // �Ҽ�ѡ�еĽڵ�

	//****** Details********
	QButtonGroup * dimensionGroupBtn;
	QRadioButton * ThreeDBtn;
	QRadioButton * TwoDBtn;
	QButtonGroup * fieldGroupBtn;
	QRadioButton * ExBtn;
	QRadioButton * EyBtn;
	QRadioButton * EzBtn;
	QRadioButton * HxBtn;
	QRadioButton * HyBtn;
	QRadioButton * HzBtn;
	QButtonGroup * pmGroupBtn;
	QRadioButton * magnitudeBtn;
	QRadioButton * phaseBtn;
	QButtonGroup * powerGroupBtn;
	QRadioButton * linearBtn;
	QRadioButton * dbBtn;

	QLabel * effLabelVal;
	QLabel * scaleLabelVal;
	QLabel * vecLabelVal;

	vector<QTreeWidgetItem*> mirrorTreeWidgetItem;

	MyData * myData;
	Mirror * tempMirror;
	Restriction * tempRestriction;

	// ����ָ��
	ApertureFieldWidget * apertureFieldWidget;
	RestrictionWidget * restrictionWidget;
	GaussianWidget * gaussianWidget;
	//ParabolicCylinderWidget * parabolicCylinderWidget;

	GraphTransWidget * tempWidget;
	FDTDProgressDialog * FDTDprogressDialog;
	PhsCorProgressDialog * phsCorprogressDialog;
	PVVAProgressDialog * PVVAprogressDialog;
	FDTDRadiator *FDTDradiator;
	Field *radiatorField; //��������������ĳ�

	bool isExistenceOpenWin;
	bool isNew;
	int fieldNum;

};

#endif // MAINWINDOW_H
