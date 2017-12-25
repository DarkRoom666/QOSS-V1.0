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

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>

#include "MyData.h"

#include "Qt/include/ParaboloidWidget.h"

using namespace userInterface;
class mainWindow : public QMainWindow
{
	Q_OBJECT

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

	void createProject();

	// 
	void isNeedSave();

	void updateVtk();


private slots:

	// �˵���Ӧ����
	void openFile();
	void newFile();
	void on_isShowBox();

	void on_modifyingMirror();
	void on_modifyParameters();

	void on_createParaboloid();
	void toReceiveParaboloid(int);

	// ------------------- ����Ҽ����� ----------------------------------
	void on_treeWidget_ContextMenuRequested(QPoint pos);// �Ҽ��˵�
	void on_treeWidget_leftPressed(QTreeWidgetItem *item, int column);

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

	// �Ҽ�
	QAction * modifyingMirrorAction;
	QAction * restrictionAction;
	QAction * modifyParametersAction;

	//------------- TreeWidgetItem------------------- 
	QTreeWidgetItem * definitionsTreeItem;
	QTreeWidgetItem * variablesTreeItem;
	QTreeWidgetItem * modelTreeItem;
	QTreeWidgetItem * geometryTreeItem;
	QTreeWidgetItem * sourceTreeItem;
	QTreeWidgetItem * lightTreeItem;
	QTreeWidgetItem * fieldTreeItem;

	QTreeWidgetItem * rightSelectItem;
	QTreeWidgetItem * leftSelectItem; // �Ҽ�ѡ�еĽڵ�

	vector<QTreeWidgetItem*> mirrorTreeWidgetItem;

	MyData * myData;
	Mirror * tempMirror;

	// ����ָ��
	ParaboloidWidget *paraboloidWidget;

	bool isExistenceOpenWin;

};

#endif // MAINWINDOW_H
