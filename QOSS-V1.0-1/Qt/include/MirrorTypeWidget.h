/*
*	created by liyun 2017/12/8
*   function ���� ѡ��������
*   y^2 = 4Fx  // fΪ����
*   R Ϊy�����ֵ
*   version 1.0
*/

#ifndef MIRRORTYPEWIDGET_H
#define MIRRORTYPEWIDGET_H

#include <QtWidgets/QDialog>
#include <QTabWidget>
#include <QPushButton>


namespace  userInterface {
	class MirrorTypeWidget : public QDialog
	{
		Q_OBJECT

	public:
		MirrorTypeWidget(QWidget *parent = 0);
		~MirrorTypeWidget();

	private slots:
		
	private:

		QTabWidget * tabWidget;

		//page1
		QPushButton *planeBtn;
		QPushButton *ellipsoidBtn;
		QPushButton *paraboloidBtn;
		QPushButton *parabolicCylinderBtn;
	};
}



#endif // MIRRORTYPEWIDGET_H
