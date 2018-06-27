#include "..\include\WaveguideWidget.h"
#include <QVBoxLayout>
#include <QFileDialog>
#include "../MyData.h"


using namespace userInterface;

WaveguideWidget::WaveguideWidget(QWidget *parent)
	: QMainWindow(parent)
{
	setCentralWidget(&widget);

	renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(1.0, 1.0, 1.0);

	auto window = widget.GetRenderWindow();
	window->AddRenderer(renderer);

	interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(window);

	auto style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	interactor->SetInteractorStyle(style);
	interactor->Initialize();

	renderer->ResetCamera();
	window->Render();

	rigthWidget = new QDockWidget;
	
	rigthWidget = new QDockWidget(this);
	rigthWidget->setFeatures(QDockWidget::NoDockWidgetFeatures);
	rigthWidget->setAllowedAreas(Qt::RightDockWidgetArea );
	addDockWidget(Qt::RightDockWidgetArea, rigthWidget);

	QWidget * qw = new QWidget(this);
	rigthWidget->setTitleBarWidget(qw);

	numBaseComboBox = new QComboBox();
	for (int i = 1; i <= 5; ++i)
	{
		numBaseComboBox->addItem(QString::number(i));
	}
	numBaseComboBox->setCurrentIndex(2);
	connect(numBaseComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(on_numBaseComboBox(int)));

	numBaseLabel = new QLabel(tr("N * N"));
	sourcePathLabel = new QLabel(tr("source path"));
	sourcePathLineEidt = new QLineEdit();
	sourcePathLineEidt->setEnabled(false);
	sourcePathBtn = new QPushButton(tr("Browse..."));
	connect(sourcePathBtn, SIGNAL(clicked()), this, SLOT(on_sourcePathBtn()));

	StartNumColumnlabel = new QLabel(tr("Start reading from Column number"));
	StartNumlabel = new QLabel(tr("Start reading from line number"));
	StartNumColumnEidt = new QLineEdit;
	StartNumLineEidt = new QLineEdit;

	StartNumLineEidt->setText(tr("1"));
	StartNumColumnEidt->setText(tr("1"));

	UNumlabel = new QLabel(tr("Number of points along U"));
	VNumlabel = new QLabel(tr("Number of points along V"));

	UNumLineEidt = new QLineEdit(tr("51"));
	VNumLineEidt = new QLineEdit(tr("51"));

	dsLabel = new QLabel(tr("source path"));
	dsLineEidt = new QLineEdit(tr("0.05"));

	datadeflabel = new QLabel(tr("Aperture data definition"));
	datadefComboBox = new QComboBox;
	datadefComboBox->addItem(tr("Amplitude & Phase"));
	datadefComboBox->addItem(tr("Real & Imaginary"));

	okBtn = new QPushButton(tr("OK"));
	connect(okBtn, SIGNAL(clicked()), this, SLOT(on_OK1Btn()));

	//QGridLayout * layoutA = new QGridLayout;
	//layoutA->addWidget(numBaseLabel, 0, 0);
	//layoutA->addWidget(numBaseComboBox, 0, 1);

	QGridLayout * layoutB = new QGridLayout;
	layoutB->addWidget(sourcePathLabel, 0, 0);
	layoutB->addWidget(sourcePathLineEidt, 0, 1);
	layoutB->addWidget(sourcePathBtn, 0, 2);

	QGridLayout * layoutC = new QGridLayout;
	layoutC->addWidget(StartNumColumnlabel, 0, 0);
	layoutC->addWidget(StartNumColumnEidt, 0, 1);
	layoutC->addWidget(StartNumlabel, 1, 0);
	layoutC->addWidget(StartNumLineEidt, 1, 1);
	layoutC->addWidget(UNumlabel, 2, 0);
	layoutC->addWidget(VNumlabel, 3, 0);
	layoutC->addWidget(UNumLineEidt, 2, 1);
	layoutC->addWidget(VNumLineEidt, 3, 1);
	layoutC->addWidget(dsLabel, 4, 0);
	layoutC->addWidget(dsLineEidt, 4, 1);
	layoutC->addWidget(datadeflabel, 5, 0);
	layoutC->addWidget(datadefComboBox, 5, 1);
	layoutC->addWidget(okBtn, 6, 1);

	QVBoxLayout * tabLayoutSource = new QVBoxLayout;
	//tabLayoutSource->addLayout(layoutA);
	tabLayoutSource->addLayout(layoutB);
	tabLayoutSource->addLayout(layoutC);

	baseGroupBox = new QGroupBox;
	baseGroupBox->setTitle(tr("Base parameter"));
	baseGroupBox->setLayout(tabLayoutSource);

	// 
	columnGapLabel = new QLabel(tr("column gap"));
	rowGapLabel = new QLabel(tr("row gap"));
	columnGapLineEdit = new QLineEdit(tr("0"));
	rowGapLineEdit = new QLineEdit(tr("0"));

	QGridLayout * layout2 = new QGridLayout;
	layout2->addWidget(numBaseLabel, 0, 0);
	layout2->addWidget(numBaseComboBox, 0, 1);
	layout2->addWidget(columnGapLabel, 1, 0);
	layout2->addWidget(columnGapLineEdit, 1, 1);
	layout2->addWidget(rowGapLabel, 2, 0);
	layout2->addWidget(rowGapLineEdit, 2, 1);

	layoutGroupBox = new QGroupBox;
	layoutGroupBox->setTitle(tr("Layout (unit: /ds)"));
	layoutGroupBox->setLayout(layout2);

	//
	columnLabel = new QLabel(tr("column"));
	rowLabel = new QLabel(tr("row"));
	columnComboBox = new QComboBox();
	rowComboBox = new QComboBox();
	amLabel = new QLabel(tr("Amplitude"));
	phsLabel = new QLabel(tr("Phase (deg)"));
	amLineEdit = new QLineEdit(tr("1"));
	phsLineEdit = new QLineEdit(tr("0"));

	QGridLayout * layout3 = new QGridLayout;
	layout3->addWidget(columnLabel, 0, 0);
	layout3->addWidget(columnComboBox, 0, 1);
	layout3->addWidget(rowLabel, 1, 0);
	layout3->addWidget(rowComboBox, 1, 1);
	layout3->addWidget(amLabel, 2, 0);
	layout3->addWidget(amLineEdit, 2, 1);
	layout3->addWidget(phsLabel, 3, 0);
	layout3->addWidget(phsLineEdit, 3, 1);

	unitGroupBox = new QGroupBox;
	unitGroupBox->setTitle(tr("Each unit"));
	unitGroupBox->setLayout(layout3);

	defaultBtn = new QPushButton(tr("default"));
	finishBtn = new QPushButton(tr("finish"));
	QGridLayout * layout4 = new QGridLayout;
	layout4->addWidget(defaultBtn, 0, 0);
	layout4->addWidget(finishBtn, 0, 1);

	QVBoxLayout * rigthlayout = new QVBoxLayout(this);
	rigthlayout->addWidget(baseGroupBox);
	rigthlayout->addWidget(layoutGroupBox);
	rigthlayout->addWidget(unitGroupBox);
	rigthlayout->addLayout(layout4);

	connect(finishBtn, SIGNAL(clicked()), this, SLOT(on_finishBtn()));
	connect(defaultBtn, SIGNAL(clicked()), this, SLOT(on_defaultBtn()));

	rWidget = new QWidget();
	rWidget->setLayout(rigthlayout);
	rigthWidget->setWidget(rWidget);
	on_numBaseComboBox(2);

	connect(columnComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(on_columnComboBox(int)));
	connect(rowComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(on_rowComboBox(int)));

	rowIndex = 0;
	colIndex = 0;

	field = new Field;
	waveguideRadiator = std::make_shared<WaveguideRadiator>();


}

WaveguideWidget::~WaveguideWidget()
{
}

void userInterface::WaveguideWidget::on_columnComboBox(int num)
{
	// save 之前的
	double colData = amLineEdit->text().toDouble();
	double rowData = phsLineEdit->text().toDouble();
	amVec[rowIndex][colIndex] = colData;
	phsVec[rowIndex][colIndex] = rowData;

	// read 现在的
	colIndex = num;
	colData = amVec[rowIndex][colIndex];
	rowData = phsVec[rowIndex][colIndex];
	amLineEdit->setText(QString::number(colData));
	phsLineEdit->setText(QString::number(rowData));
}

void userInterface::WaveguideWidget::on_rowComboBox(int num)
{
	// save 之前的
	double colData = amLineEdit->text().toDouble();
	double rowData = phsLineEdit->text().toDouble();
	amVec[rowIndex][colIndex] = colData;
	phsVec[rowIndex][colIndex] = rowData;

	// read 现在的
	rowIndex = num;
	colData = amVec[rowIndex][colIndex];
	rowData = phsVec[rowIndex][colIndex];
	amLineEdit->setText(QString::number(colData));
	phsLineEdit->setText(QString::number(rowData));
}

void userInterface::WaveguideWidget::on_sourcePathBtn()
{
	QString filename = QFileDialog::getOpenFileName(this,
		tr("Open the file"),
		"",
		tr("*.txt"));
	if (!filename.isEmpty())
	{
		sourcePathLineEidt->setText(filename);
	}
}

void userInterface::WaveguideWidget::on_OK1Btn()
{
	int N_width = UNumLineEidt->text().toInt();
	int	M_depth = VNumLineEidt->text().toInt();
	std::string fileAddress = sourcePathLineEidt->text().toStdString();
	double ds = dsLineEidt->text().toDouble();
	bool isAmPhs = false;
	if (datadefComboBox->currentIndex() == 0)
		isAmPhs = true;
	int beginY = StartNumLineEidt->text().toInt();
	int beginX = StartNumColumnEidt->text().toInt();

	waveguideRadiator->setBasePara(fileAddress, N_width, M_depth, beginY-1, beginX-1, ds, isAmPhs);
	waveguideRadiator->readData();
	waveguideRadiator->setEachUnit(amVec, phsVec);
	waveguideRadiator->genField(field);
	field->updateData();
	renderer->AddActor(field->getActor());

	auto window = widget.GetRenderWindow();
	window->Render();
	
}

void userInterface::WaveguideWidget::on_finishBtn()
{
	double colData = amLineEdit->text().toDouble();
	double rowData = phsLineEdit->text().toDouble();
	amVec[rowIndex][colIndex] = colData;
	phsVec[rowIndex][colIndex] = rowData;

	int gapCol = columnGapLineEdit->text().toInt();
	int gapRow = rowGapLineEdit->text().toInt();
	int numNN = numBaseComboBox->currentIndex() + 1;
	waveguideRadiator->setLayout(gapCol, gapRow, numNN);
	waveguideRadiator->setEachUnit(amVec, phsVec);
	waveguideRadiator->genField(field);
	field->updateData();
	//renderer->AddActor(field->getActor());

	auto window = widget.GetRenderWindow();
	window->Render();
}

void userInterface::WaveguideWidget::on_defaultBtn()
{
	numBaseComboBox->setCurrentIndex(2); // 会自动调用 on_numBaseComboBox
	columnGapLineEdit->setText(tr("0"));
	rowGapLineEdit->setText(tr("0"));
	on_finishBtn();
}

void WaveguideWidget::on_numBaseComboBox(int num)
{
	disconnect(columnComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(on_columnComboBox(int)));
	disconnect(rowComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(on_rowComboBox(int)));

	columnComboBox->clear();
	rowComboBox->clear();
	for (int i = 0; i <= num; ++i)
	{
		columnComboBox->addItem(QString::number(i + 1));
		rowComboBox->addItem(QString::number(i + 1));
	}
	amVec.clear();
	phsVec.clear();
	amVec.resize(num + 1);
	phsVec.resize(num + 1);
	for (int i = 0; i <= num; ++i)
	{
		amVec[i].resize(num + 1);
		phsVec[i].resize(num + 1);
	}
	for (int i = 0; i <= num; ++i)
		for (int j = 0; j <= num; ++j)
		{
			amVec[i][j] = 1;
			phsVec[i][j] = 0;
		}
	amLineEdit->setText(QString::number(1));
	phsLineEdit->setText(QString::number(0));

	connect(columnComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(on_columnComboBox(int)));
	connect(rowComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(on_rowComboBox(int)));

	rowIndex = 0;
	colIndex = 0;
}



