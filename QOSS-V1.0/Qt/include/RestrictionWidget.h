/*
*	created by liyun 2017/12/26  
*   function ������������
*   version 1.0
*/
#ifndef RESTRICTIONWIDGET_H
#define RESTRICTIONWIDGET_H

#include <QtWidgets/QDialog>
#include <QTabWidget>
#include <QGroupBox>
#include <QGridLayout>
#include <Qlabel>
#include <QLineEdit>
#include <QPushButton>
#include <QComboBox>

#include "Qt/include/GraphTransWidget.h"
#include "../VTK/include/Restriction.h"


namespace  userInterface {
	class RestrictionWidget : public GraphTransWidget
	{
		Q_OBJECT

	public:
		RestrictionWidget(QWidget *parent = 0);
		~RestrictionWidget();
		void setRestriction(Restriction*);

	private slots:
		void on_radiusChange(QString);
		void on_focusChange(QString);
		void on_typeComboBox(int index);

	private:
		//page1
		QTabWidget * tabWidget;
		QGroupBox * defGroupBox;
		QGroupBox * baseGroupBox;
		QGroupBox * dimGroupBox;

		QLabel * radiuslabel;
		QLabel * depthlabel;

		QLineEdit * radiusLineEidt;
		QLineEdit * depthLineEidt;

		QLabel * label;
		QLineEdit * nameLineEidt;

		QLabel * typeLabel;
		QComboBox * typeComboBox;

		Restriction * restriction;
	};
}


#endif // RESTRICTIONWIDGET_H
