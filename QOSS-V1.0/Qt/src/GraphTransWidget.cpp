#include "Qt/include/GraphTransWidget.h"

#include <QMessageBox>
using namespace userInterface;

GraphTransWidget::GraphTransWidget(QWidget *parent, int wayButton)
	: QDialog(parent)
{
	this->wayButton = wayButton;
}

GraphTransWidget::~GraphTransWidget()
{

}

void GraphTransWidget::addDefGroupBox(QGroupBox * _defGroupBox, QString filename)
{
	// page1
	//defGroupBox
	QLabel * imgLlabel;
	imgLlabel = new QLabel;
	imgLlabel->setPixmap(QPixmap(filename));

	QGridLayout * layout1 = new QGridLayout;
	layout1->addWidget(imgLlabel);

	//defGroupBox = new QGroupBox;
	_defGroupBox->setLayout(layout1);
	_defGroupBox->setTitle(tr("Definition methods"));
}

void GraphTransWidget::addBaseGroupBox(QGroupBox * _baseGroupBox)
{
	//baseGroupBox
	Ulabel = new QLabel(tr("U"));
	Vlabel = new QLabel(tr("V"));
	Nlabel = new QLabel(tr("N"));

	ULineEidt = new QLineEdit;
	VLineEidt = new QLineEdit;
	NLineEidt = new QLineEdit;

	ULineEidt->setText(tr("0.0"));
	VLineEidt->setText(tr("0.0"));
	NLineEidt->setText(tr("0.0"));

	connect(ULineEidt, SIGNAL(textChanged(QString)),
		this, SLOT(on_GraphChange(QString)));
	connect(VLineEidt, SIGNAL(textChanged(QString)),
		this, SLOT(on_GraphChange(QString)));
	connect(NLineEidt, SIGNAL(textChanged(QString)),
		this, SLOT(on_GraphChange(QString)));

	QGridLayout * layout2 = new QGridLayout;
	layout2->addWidget(Ulabel, 0, 0);
	layout2->addWidget(Vlabel, 1, 0);
	layout2->addWidget(Nlabel, 2, 0);
	layout2->addWidget(ULineEidt, 0, 1);
	layout2->addWidget(VLineEidt, 1, 1);
	layout2->addWidget(NLineEidt, 2, 1);

	_baseGroupBox->setTitle(tr("Definition methods"));
	_baseGroupBox->setLayout(layout2);
}

void GraphTransWidget::addRotateWidget(QWidget * RotateWidget, QString filename)
{
	// page2
	//imgGroupBox
	imgLlabel1 = new QLabel;
	imgLlabel1->setPixmap(QPixmap(filename));

	QGridLayout * layout5 = new QGridLayout;
	layout5->addWidget(imgLlabel1);

	imgGroupBox = new QGroupBox;
	imgGroupBox->setLayout(layout5);
	imgGroupBox->setTitle(tr("Definition methods"));

	//RotateGroupBox
	xlabel = new QLabel(tr("X"));
	ylabel = new QLabel(tr("Y"));
	zlabel = new QLabel(tr("Z"));

	xRotateLineEidt = new QLineEdit;
	yRotateLineEidt = new QLineEdit;
	zRotateLineEidt = new QLineEdit;
	connect(xRotateLineEidt, SIGNAL(textChanged(QString)),
		this, SLOT(on_GraphChange(QString)));
	connect(yRotateLineEidt, SIGNAL(textChanged(QString)),
		this, SLOT(on_GraphChange(QString)));
	connect(zRotateLineEidt, SIGNAL(textChanged(QString)),
		this, SLOT(on_GraphChange(QString)));

	xRotateLineEidt->setText(tr("1.0"));
	yRotateLineEidt->setText(tr("0.0"));
	zRotateLineEidt->setText(tr("0.0"));

	QGridLayout * layout6 = new QGridLayout;
	layout6->addWidget(xlabel, 0, 0);
	layout6->addWidget(ylabel, 1, 0);
	layout6->addWidget(zlabel, 2, 0);
	layout6->addWidget(xRotateLineEidt, 0, 1);
	layout6->addWidget(yRotateLineEidt, 1, 1);
	layout6->addWidget(zRotateLineEidt, 2, 1);

	RotateGroupBox = new QGroupBox;
	RotateGroupBox->setTitle(tr("Axis of rotation"));
	RotateGroupBox->setLayout(layout6);

	//thetaGroupBox
	thetalabel = new QLabel(tr("degree"));
	thetaLineEidt = new QLineEdit;
	thetaLineEidt->setText(tr("0.0"));

	connect(thetaLineEidt, SIGNAL(textChanged(QString)),
		this, SLOT(on_GraphChange(QString)));

	QGridLayout * layout7 = new QGridLayout;
	layout7->addWidget(thetalabel, 0, 0);
	layout7->addWidget(thetaLineEidt, 0, 1);

	thetaGroupBox = new QGroupBox;
	thetaGroupBox->setTitle(tr("Rotation angle"));
	thetaGroupBox->setLayout(layout7);

	QVBoxLayout * tabLayout2;
	tabLayout2 = new QVBoxLayout;
	tabLayout2->addWidget(imgGroupBox);
	tabLayout2->addWidget(RotateGroupBox);
	tabLayout2->addWidget(thetaGroupBox);

	RotateWidget->setLayout(tabLayout2);
}

void GraphTransWidget::addBtn(QGridLayout * _layoutbt, int wayButton)
{
	this->wayButton = wayButton;
	if (wayButton == 1)
	{
		createbtn = new QPushButton(tr("Ok"));
		addbtn = new QPushButton(tr("Apply"));
		closebtn = new QPushButton(tr("Cancle"));

		connect(createbtn, SIGNAL(clicked()), this, SLOT(buttonOk()));
		connect(addbtn, SIGNAL(clicked()), this, SLOT(buttonApply()));
		connect(closebtn, SIGNAL(clicked()), this, SLOT(buttonClose()));
	}
	else
	{
		createbtn = new QPushButton(tr("Create"));
		addbtn = new QPushButton(tr("Add"));
		closebtn = new QPushButton(tr("Close"));

		connect(createbtn, SIGNAL(clicked()), this, SLOT(buttonCreate()));
		connect(addbtn, SIGNAL(clicked()), this, SLOT(buttonAdd()));
		connect(closebtn, SIGNAL(clicked()), this, SLOT(buttonClose()));
	}

	_layoutbt->addWidget(createbtn, 0, 0);
	_layoutbt->addWidget(addbtn, 0, 1);
	_layoutbt->addWidget(closebtn, 0, 2);
}

void GraphTransWidget::buttonApply()
{
	emit sendData(6);
	
}

void GraphTransWidget::buttonOk()
{
	emit sendData(3);
}

