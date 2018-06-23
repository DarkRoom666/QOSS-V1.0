/*
*	created by liyun 2018/6/23
*   function ´´½¨Íø¸ñ
*   version 1.0
*/

#ifndef MESHWIDGET_H
#define MESHWIDGET_H

#include <QtWidgets/QDialog>
#include <QComboBox>
#include <QLabel>
#include <QPushButton>
#include <QProgressBar>
#include <QGroupBox>


namespace  userInterface {
	class MeshWidget : public QDialog
	{
		Q_OBJECT

	public:
		MeshWidget(QWidget *parent = 0);
		~MeshWidget();



	private slots:
		void on_Start();
		void on_stop();
		void setMainValue(int);

	protected:
		void closeEvent(QCloseEvent *event);
	private:

		QGroupBox * defGroupBox;
		QComboBox * accuracyComboBox;
		QLabel* accuracyLabel;
		QLabel* meshLabel;
		QLabel* meshNumLabel;

		QPushButton *startBtn;

		QProgressBar *mainBar;
		QProgressBar *slaverBar;
		QLabel * txtLabel;

		QPushButton *closeBtn;

		int MirrorNum;

	};
}



#endif // MESHWIDGET_H
