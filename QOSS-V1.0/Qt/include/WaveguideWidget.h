/*
*	created by liyun 2018/6/23
*   function 创建网格
*   version 1.0
*/

#ifndef WAVEGUIDEWIDGET_H
#define WAVEGUIDEWIDGET_H

#include <QtWidgets/QMainWindow>
#include <QComboBox>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QGroupBox>
#include <QVTKWidget.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <QDockWidget>
#include <vector>
#include <memory>

#include "../VTK/include/WaveguideRadiator.h"
#include "../VTK/include/Field.h"


namespace  userInterface {
	class WaveguideWidget : public QMainWindow
	{
		Q_OBJECT

	public:
		WaveguideWidget(QWidget *parent = 0);
		~WaveguideWidget();


	private slots:
		void on_numBaseComboBox(int);
		void on_columnComboBox(int);
		void on_rowComboBox(int);
		void on_sourcePathBtn();
		void on_OK1Btn();


	private:
		QVTKWidget widget; // vtk 显示窗口
		vtkSmartPointer<vtkRenderWindowInteractor> interactor; // 交互
		vtkSmartPointer<vtkRenderer> renderer;

		QDockWidget * rigthWidget;
		QWidget * rWidget;

		QGroupBox * baseGroupBox;
		QComboBox * numBaseComboBox;
		QLabel* numBaseLabel;
		QLabel* sourcePathLabel;
		QPushButton *sourcePathBtn;
		QLineEdit * sourcePathLineEidt;
		QLabel* dsLabel;
		QLineEdit* dsLineEidt;

		QLabel* datadeflabel;
		QComboBox * datadefComboBox;

		QGroupBox * layoutGroupBox;
		QLabel *columnGapLabel;
		QLabel *rowGapLabel;
		QLineEdit *columnGapLineEdit;
		QLineEdit *rowGapLineEdit;

		QLabel *amLabel;
		QLabel *phsLabel;
		QLineEdit *amLineEdit;
		QLineEdit *phsLineEdit;

		QGroupBox * unitGroupBox;
		QLabel *columnLabel;
		QLabel *rowLabel;
		QComboBox * columnComboBox;
		QComboBox * rowComboBox;

		QPushButton * okBtn;

		//
		QLabel *StartNumColumnlabel;
		QLabel *StartNumlabel;

		QLineEdit * StartNumLineEidt;
		QLineEdit * StartNumColumnEidt;

		QGroupBox * destinationGroupBox;
		QLabel * UNumlabel;
		QLabel * VNumlabel;

		QLineEdit * UNumLineEidt;
		QLineEdit * VNumLineEidt;

		int MirrorNum;

		std::vector<std::vector<double>> amVec;
		std::vector<std::vector<double>> phsVec;

		int rowIndex;
		int colIndex;

		Field * field;
		shared_ptr<WaveguideRadiator> waveguideRadiator;

	};
}



#endif // MESHWIDGET_H
