#include "..\include\MeshWidget.h"
#include <QVBoxLayout>
#include "../MyData.h"
#include "../include/MeshThread.h"
using namespace userInterface;

MeshWidget::MeshWidget(QWidget *parent)
	: QDialog(parent)
{
	setWindowTitle(tr("Mesh"));

	MirrorNum = MyData::getInstance()->getNumOfMirrors();
	
	accuracyLabel = new QLabel(tr("Accuracy "));
	accuracyComboBox = new QComboBox();
	accuracyComboBox->addItem(tr("Crude"));
	accuracyComboBox->addItem(tr("Standard"));
	accuracyComboBox->addItem(tr("Fine"));
	accuracyComboBox->setCurrentIndex(1);

	//meshLabel = new QLabel(tr("About the number of grids: "));
	//meshNumLabel = new QLabel;

	startBtn = new QPushButton(tr("Start"));
	connect(startBtn, SIGNAL(clicked()), this, SLOT(on_Start()));

	QGridLayout * layout = new QGridLayout;
	layout->addWidget(accuracyLabel, 0, 0);
	layout->addWidget(accuracyComboBox, 0, 1);
	layout->addWidget(startBtn, 0, 2);


	defGroupBox = new QGroupBox;
	defGroupBox->setTitle(tr("Dimensions"));
	defGroupBox->setLayout(layout);

	txtLabel = new QLabel();
	mainBar = new QProgressBar();
	mainBar->setRange(0, MirrorNum +1);
	mainBar->setMinimumWidth(500);
	slaverBar = new QProgressBar();
	slaverBar->setRange(0, 100);
	slaverBar->setMinimumWidth(500);

	QVBoxLayout * mainlayout = new QVBoxLayout(this);
	mainlayout->addWidget(defGroupBox);
	mainlayout->addWidget(txtLabel);
	mainlayout->addWidget(mainBar);
	mainlayout->addWidget(slaverBar);

	txtLabel->setHidden(true);
	mainBar->setHidden(true);
	slaverBar->setHidden(true);
}


MeshWidget::~MeshWidget()
{
}

void userInterface::MeshWidget::setMainValue(int val)
{
	mainBar->setValue(val);
	if (0 == val)
		txtLabel->setText(tr("Initializing..."));
	else if (val <= MirrorNum)
		txtLabel->setText(tr("Mesh computing Mirror") + QString::number(val + 1) + tr(" ..."));
	else
	{
		txtLabel->setText(tr("All computed, Loading Fields..."));
	}
}

void userInterface::MeshWidget::on_Start()
{
	txtLabel->setHidden(false);
	mainBar->setHidden(false);
	//slaverBar->setHidden(false);
	accuracyComboBox->setEnabled(false);
	startBtn->setEnabled(false);
	int type = accuracyComboBox->currentIndex();
	MeshThread *calThr = new MeshThread(type);
	connect(calThr, SIGNAL(sendMainValue(int)), this, SLOT(setMainValue(int)));
	//connect(calThr, SIGNAL(sendStop()), calThr, SLOT(killFDTD()));
	calThr->start();
}

void userInterface::MeshWidget::on_stop()
{
}

void userInterface::MeshWidget::closeEvent(QCloseEvent * event)
{
}


