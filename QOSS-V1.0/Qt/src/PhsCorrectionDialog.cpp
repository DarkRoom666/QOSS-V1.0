#include "..\include\PhsCorrectionDialog.h"
#include <QMessageBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
using namespace userInterface;

PhsCorrectionDialog::PhsCorrectionDialog(QWidget *parent)
	: QDialog(parent)
{
	setWindowTitle(tr("Phase Correction"));
	
	lengthLabel = new QLabel(tr("Length:"));
	lengthLineEdit = new QLineEdit(tr("0.5"));
	dsLabel = new QLabel(tr("Accuracy:"));
	dsComboBox = new QComboBox;
	dsComboBox->addItem(tr("Crude"));
	dsComboBox->addItem(tr("Standard"));
	dsComboBox->addItem(tr("Fine"));
	dsComboBox->setCurrentIndex(1);

	QGridLayout *gridLayout = new QGridLayout;
	gridLayout->addWidget(lengthLabel, 0, 0);
	gridLayout->addWidget(lengthLineEdit, 0, 1);
	gridLayout->addWidget(dsLabel, 1, 0);
	gridLayout->addWidget(dsComboBox, 1, 1);

	OKBtn = new QPushButton(tr("OK"));
	connect(OKBtn, SIGNAL(clicked()), this, SLOT(on_OK()));
	cancelBtn = new QPushButton(tr("Cancel"));
	connect(cancelBtn, SIGNAL(clicked()), this, SLOT(on_cancel()));

	QHBoxLayout *hlayout = new QHBoxLayout;
	hlayout->addSpacing(100);
	hlayout->addWidget(OKBtn);
	hlayout->addWidget(cancelBtn);

	QVBoxLayout *layout = new QVBoxLayout(this);

	layout->addLayout(gridLayout);
	layout->addLayout(hlayout);

}

PhsCorrectionDialog::~PhsCorrectionDialog()
{
}

void PhsCorrectionDialog::on_OK()
{
	accept();
}

void PhsCorrectionDialog::on_cancel()
{
	close();
}
