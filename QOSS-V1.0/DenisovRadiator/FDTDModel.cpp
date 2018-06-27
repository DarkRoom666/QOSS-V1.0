#include "../DenisovRadiator/FDTDModel.h"
#include "../DenisovRadiator/CodeJin/FieldLJ.h"

FDTDModel::FDTDModel(QWidget *parent)
	: QWidget(parent)
{
	excisopen = false;
	modelisopen = false;
	fdtdcalculated = false;
	Field3D = false;
	FDTD = new QTFDTDThread;


	FDTDModel::CreateButtons();
	QHBoxLayout *Layout = new QHBoxLayout(this);
	Layout->addWidget(&widget);
	Layout->addWidget(Buttons);
	//示例数据
	Nx = 100;
	Ny = 100;
	Nz = 100;
	dx = 0.01;
	dy = 0.01;
	dz = 0.01;
	Eps.resize(Nx*Ny*Nz);
	vector<double> xx;	xx.resize(Nx);
	vector<double> yy;	yy.resize(Ny);
	for (int i = 0; i < Nx; i++) {
		xx[i] = (i + 0.5 - Nx / 2)*dx;
	}
	for (int j = 0; j < Ny; j++) {
		yy[j] = (j + 0.5 - Ny / 2)*dy;
	}
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++){
			for (int k = 0; k < Nz; k++){
				//Eps[i + j*Nx + k*Nx*Ny] = 0.5;
				
				Eps[i + j*Nx + k*Nx*Ny] = 0.0;
				if (xx[i] * xx[i] + yy[j] * yy[j] < 0.4*0.4 && xx[i] * xx[i] + yy[j] * yy[j] > 0.3*0.3) {
					Eps[i + j*Nx + k*Nx*Ny] = 1.0;
				}
			}
		}
	}

	//越小越透明
	FDTDModel::transparency = 1;
	FDTDModel::min_gray = 0;
	FDTDModel::max_gray = 200;

	double axesScale = Nx*dx*0.5;
	resize(1400, 800);
	//设置三维数组
	auto img = vtkSmartPointer<vtkImageData>::New();
	auto info = vtkSmartPointer<vtkInformation>::New();
	//auto info3 = img->GetInformation();//此处也可以从img里面获取information，如
	//果information为空，则会创建一个新的，相当于调用New之后再SetInformation
	img->SetDimensions(Nx, Ny, Nz);
	img->SetSpacing(dx, dy, dz);
	img->SetOrigin(0, 0, 0);
	img->SetScalarType(VTK_UNSIGNED_CHAR, info);
	img->SetNumberOfScalarComponents(1, info);
	img->AllocateScalars(info);
	unsigned char *ptr = (unsigned char *)img->GetScalarPointer();
	for (int i = 0; i < Nx*Ny*Nz; i++){
		*ptr++ = float(Eps[i]) * max_gray;
	}

	volumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	volumeMapper->SetInputData(img);

	vtkSmartPointer<vtkVolumeProperty> volumeProperty
		= vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->SetInterpolationTypeToLinear();
	volumeProperty->SetAmbient(0.1);
	volumeProperty->SetDiffuse(0.1);
	volumeProperty->SetSpecular(0.1);

	vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity =
		vtkSmartPointer<vtkPiecewiseFunction>::New();
	compositeOpacity->ClampingOn();
	compositeOpacity->AddPoint(min_gray, 0);
	compositeOpacity->AddPoint(max_gray, transparency);

	volumeProperty->SetScalarOpacity(compositeOpacity); //设置不透明度传输函数

	vtkSmartPointer<vtkColorTransferFunction> color =
		vtkSmartPointer<vtkColorTransferFunction>::New();

	float temp = (255 - min_gray) / 7;
	color->AddRGBPoint(min_gray, 0.0, 0.0, 0.0);
	color->AddRGBPoint(min_gray + temp, 1.0, 0.0, 0.0);
	color->AddRGBPoint(min_gray + temp * 2, 1.0, 1.0, 0.0);
	color->AddRGBPoint(min_gray + temp * 3, 0.0, 1.0, 0.0);
	color->AddRGBPoint(min_gray + temp * 4, 0.0, 1.0, 1.0);
	color->AddRGBPoint(min_gray + temp * 5, 0.0, 0.0, 1.0);
	color->AddRGBPoint(min_gray + temp * 6, 1.0, 0.0, 1.0);
	color->AddRGBPoint(min_gray + temp * 7, 1.0, 1.0, 1.0);
	volumeProperty->SetColor(color);

	volumeProperty->SetColor(color);

	volume = vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(volumeMapper);
	volume->SetProperty(volumeProperty);

	renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(0.0, 0.0, 0.0);
	renderer->AddVolume(volume);
	

	/*
	//视角
	vtkCamera *aCamera = vtkCamera::New();
	aCamera->SetViewUp(0, 0, 1);//视角
	aCamera->SetPosition(0, -3 * axesScale, 0);
	aCamera->SetFocalPoint(0, 0, 0);
	aCamera->ComputeViewPlaneNormal();
	
	renderer->SetActiveCamera(aCamera);
	renderer->ResetCamera();
	*/
	auto window = widget.GetRenderWindow();
	window->AddRenderer(renderer);
	window->Render();

	//交互
	interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(window);

	auto style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	interactor->SetInteractorStyle(style);
	interactor->Initialize();

	//坐标轴
	axes = vtkSmartPointer<vtkAxesActor>::New();
	axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(1, 0, 0);
	axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0, 2, 0);
	axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0, 0, 3);
	axes->SetConeRadius(0.3);	axes->SetConeResolution(20);

	widget1 = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
	widget1->SetOutlineColor(0.9300, 0.5700, 0.1300);
	widget1->SetOrientationMarker(axes);
	widget1->SetInteractor(interactor);
	widget1->SetViewport(0.0, 0.0, 0.25, 0.25);
	widget1->SetEnabled(1);
	widget1->InteractiveOff();

	//Slots
	connect(OpenModel, SIGNAL(clicked()), this, SLOT(OpenModelFile()));
	connect(OpenExc, SIGNAL(clicked()), this, SLOT(OpenExcFile()));

	//Slots
	connect(SetFDTD, SIGNAL(clicked()), this, SLOT(SetFDTDprocess()));
	connect(CalculateFDTD, SIGNAL(clicked()), this, SLOT(RunFDTDprocess()));
	//
	connect(FDTDModel::FDTD, SIGNAL(sendMainValue(int)), this, SLOT(recieveInt(int)));
	connect(FDTDModel::FDTD, SIGNAL(sendSlaverValue(float)), this, SLOT(recieveFloat(float)));
}

FDTDModel::~FDTDModel() {

}

void FDTDModel::CreateButtons(void) {

	//默认，命名 按钮
	OpenModel = new QPushButton(tr("OpenModel"));	OpenModel->setMaximumWidth(150);
	OpenExc = new QPushButton(tr("OpenExc"));		OpenExc->setMaximumWidth(150);
	SetFreList = new QPushButton(tr("SetFreList"));	SetFreList->setMaximumWidth(150);
	SetFreList->setEnabled(false);
	SetFDTD = new QPushButton(tr("SetFDTD"));		SetFDTD->setMaximumWidth(150);
	SetFDTD->setEnabled(false);
	CalculateFDTD = new QPushButton(tr("Run FDTD")); CalculateFDTD->setMaximumWidth(150);
	CalculateFDTD->setEnabled(false);
	SetAperture = new QPushButton(tr("SetAperture"));	SetAperture->setMaximumWidth(150);
	SetAperture->setEnabled(false);
	CalAperture = new QPushButton(tr("CalAperture"));	CalAperture->setMaximumWidth(150);
	CalAperture->setEnabled(false);
	//参数列表
	FreqLabel = new QLabel(tr("Freq(GHz)"));	FreqLabel->setMaximumWidth(150);
	FreqEdit = new QLineEdit(tr("140"));		FreqEdit->setMaximumWidth(150);
	FreqEdit->setEnabled(false);
	
	dxLabel = new QLabel(tr("dx in m"));		dxLabel->setMaximumWidth(150);
	dxEdit = new QLineEdit(tr("dx"));			dxEdit->setMaximumWidth(150);
	dxEdit->setEnabled(false);
	
	dyLabel = new QLabel(tr("dy in m"));		dyLabel->setMaximumWidth(150);
	dyEdit = new QLineEdit(tr("dy"));			dyEdit->setMaximumWidth(150);
	dyEdit->setEnabled(false);
	
	dzLabel = new QLabel(tr("dz in m"));		dzLabel->setMaximumWidth(150);
	dzEdit = new QLineEdit(tr("dz"));			dzEdit->setMaximumWidth(150);
	dzEdit->setEnabled(false);

	NxLabel = new QLabel(tr("Nx"));				NxLabel->setMaximumWidth(150);
	NxEdit = new QLineEdit("Nx");				NxEdit->setMaximumWidth(150);
	NxEdit->setEnabled(false);

	NyLabel = new QLabel(tr("Ny"));				NyLabel->setMaximumWidth(150);
	NyEdit = new QLineEdit("Ny");				NyEdit->setMaximumWidth(150);
	NyEdit->setEnabled(false);

	NzLabel = new QLabel(tr("Nz"));				NzLabel->setMaximumWidth(150);
	NzEdit = new QLineEdit("Nz");				NzEdit->setMaximumWidth(150);
	NzEdit->setEnabled(false);

	OMPLabel = new QLabel(tr("OMP Num"));		OMPLabel->setMaximumWidth(150);
	OMPEdit = new QLineEdit(tr("4"));			OMPEdit->setMaximumWidth(150);
	//进度
	Message = new QLabel(tr("Progress Message"));
	progressBar = new QProgressBar;
	progressBar->setMaximumWidth(300);
	progressBar->setRange(0, 100);
	progressBar->setValue(0);
	progressBar->setHidden(true);

	Buttons = new QGroupBox(tr("Operations")); Buttons->setMaximumWidth(320);
	QVBoxLayout *ButtonV = new QVBoxLayout(Buttons);

	//顶部区域
	QHBoxLayout *H1 = new QHBoxLayout;
	H1->addWidget(OpenModel);	H1->addWidget(OpenExc);
	ButtonV->addLayout(H1);
	//模型显示透明调节区
	transparencyLabel = new QLabel(tr("transparency:"));
	transparencyLineEdit = new QLineEdit(tr("1"));
	QGridLayout * transparencyLabelLayout = new QGridLayout;
	transparencyLabelLayout->addWidget(transparencyLabel, 0, 0);
	transparencyLabelLayout->addWidget(transparencyLineEdit, 0, 1);
	transparencyScrollBar = new QScrollBar(Qt::Horizontal);
	transparencyScrollBar->setRange(0, 100);
	transparencyScrollBar->setSingleStep(1);
	transparencyScrollBar->setValue(100);
	connect(transparencyScrollBar, SIGNAL(valueChanged(int)),
		this, SLOT(transparencyChanged(int)));
	QVBoxLayout * layoutTrans = new QVBoxLayout;
	layoutTrans->addLayout(transparencyLabelLayout);
	layoutTrans->addWidget(transparencyScrollBar);
	ButtonV->addLayout(layoutTrans);
	//参数区
	QVBoxLayout *V1 = new QVBoxLayout;
	QVBoxLayout *V2 = new QVBoxLayout;
	V1->addWidget(FreqLabel);	V2->addWidget(FreqEdit);
	V1->addWidget(dxLabel);		V2->addWidget(dxEdit);
	V1->addWidget(dyLabel);		V2->addWidget(dyEdit);
	V1->addWidget(dzLabel);		V2->addWidget(dzEdit);
	V1->addWidget(NxLabel);		V2->addWidget(NxEdit);
	V1->addWidget(NyLabel);		V2->addWidget(NyEdit);
	V1->addWidget(NzLabel);		V2->addWidget(NzEdit);
	V1->addWidget(OMPLabel);	V2->addWidget(OMPEdit);
	QHBoxLayout *H2 = new QHBoxLayout;
	H2->addLayout(V1);	H2->addLayout(V2);
	ButtonV->addLayout(H2);
	//设置计算区
	ButtonV->addWidget(SetFreList);
	QHBoxLayout *H3 = new QHBoxLayout;
	H3->addWidget(SetFDTD);	H3->addWidget(CalculateFDTD);
	ButtonV->addLayout(H3);
	QHBoxLayout *H4 = new QHBoxLayout;
	H4->addWidget(SetAperture);	H4->addWidget(CalAperture);
	ButtonV->addLayout(H4);
	//进度区
	QVBoxLayout *V3 = new QVBoxLayout;
	V3->addWidget(Message);
	V3->addWidget(progressBar);
	ButtonV->addLayout(V3);

}

void FDTDModel::OpenModelFile(void) {
	QString filename = QFileDialog::getOpenFileName(this,
		tr("Open the file"),
		"",
		tr("*.dat"));
	if (!filename.isEmpty())
	{
		//读文件
		FILE* modelin;
		modelfile = filename.toStdString();
		modelin = fopen(modelfile.c_str(), "rb");
		int readi; float readf;
		fread(&readi, sizeof(int), 1, modelin);	FDTDModel::Nx = readi;
		fread(&readi, sizeof(int), 1, modelin);	FDTDModel::Ny = readi;
		fread(&readi, sizeof(int), 1, modelin);	FDTDModel::Nz = readi;
		fread(&readi, sizeof(int), 1, modelin); FDTDModel::Shiftx = readi;
		fread(&readi, sizeof(int), 1, modelin);	FDTDModel::Shifty = readi;
		
		fread(&readf, sizeof(float), 1, modelin); FDTDModel::dx = readf;
		fread(&readf, sizeof(float), 1, modelin); FDTDModel::dy = readf;
		fread(&readf, sizeof(float), 1, modelin); FDTDModel::dz = readf;

		fread(&readf, sizeof(float), 1, modelin); FDTDModel::cx = readf;
		fread(&readf, sizeof(float), 1, modelin); FDTDModel::cy = readf;
		fread(&readf, sizeof(float), 1, modelin); FDTDModel::cz = readf;

		//改变参数显示
		FDTDModel::NxEdit->setText(QString::number(Nx, 10, 0));
		FDTDModel::NyEdit->setText(QString::number(Ny, 10, 0));
		FDTDModel::NzEdit->setText(QString::number(Nz, 10, 0));
		FDTDModel::dxEdit->setText(QString::number(dx, 10, 6));
		FDTDModel::dyEdit->setText(QString::number(dy, 10, 6));
		FDTDModel::dzEdit->setText(QString::number(dz, 10, 6));
		//设置判断
		modelisopen = true;
		if (modelisopen&excisopen)
			FDTDModel::SetFDTD->setEnabled(true);

		Eps.resize(FDTDModel::Nx*FDTDModel::Ny*FDTDModel::Nz);
		for (int k = 0; k < FDTDModel::Nz; k++) {
			for (int j = 0; j < FDTDModel::Ny; j++) {
				for (int i = 0; i <FDTDModel::Nx; i++) {
					
					fread(&readf, sizeof(float), 1, modelin);
					Eps[i + j*Nx + k*Nx*Ny] = readf;
				}
			}
		}
		fclose(modelin);

		auto img = vtkSmartPointer<vtkImageData>::New();
		auto info = vtkSmartPointer<vtkInformation>::New();
		//auto info3 = img->GetInformation();//此处也可以从img里面获取information，如
		//果information为空，则会创建一个新的，相当于调用New之后再SetInformation
		img->SetDimensions(FDTDModel::Nx, FDTDModel::Ny, FDTDModel::Nz);
		img->SetSpacing(FDTDModel::dx, FDTDModel::dy, FDTDModel::dz);
		//注意是Origin哦~
		img->SetOrigin(FDTDModel::cx-dx*Nx*0.5, FDTDModel::cy-dy*Ny*0.5, FDTDModel::cz-dz*Nz*0.5);
		img->SetScalarType(VTK_UNSIGNED_CHAR, info);
		img->SetNumberOfScalarComponents(1, info);
		img->AllocateScalars(info);

		unsigned char *ptr = (unsigned char *)img->GetScalarPointer();
		for (int i = 0; i < Nx*Ny*Nz; i++) {
			*ptr++ = float(Eps[i]) *max_gray;
		}

		//volumeMapper = vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
		volumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
		volumeMapper->SetInputData(img);

		vtkSmartPointer<vtkVolumeProperty> volumeProperty
			= vtkSmartPointer<vtkVolumeProperty>::New();
		volumeProperty->SetInterpolationTypeToLinear();
		volumeProperty->ShadeOn();  //打开或者关闭阴影测试
		volumeProperty->SetAmbient(0.1);
		volumeProperty->SetDiffuse(0.1);
		volumeProperty->SetSpecular(0.1);

		vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity =
			vtkSmartPointer<vtkPiecewiseFunction>::New();
		compositeOpacity->ClampingOn();
		compositeOpacity->AddPoint(min_gray, 0);
		compositeOpacity->AddPoint(max_gray, transparency);

		volumeProperty->SetScalarOpacity(compositeOpacity); //设置不透明度传输函数


		vtkSmartPointer<vtkColorTransferFunction> color =
			vtkSmartPointer<vtkColorTransferFunction>::New();
		/*
		float temp = (255 - min_gray) / 7;
		color->AddRGBPoint(min_gray, 0.0, 0.0, 0.0);
		color->AddRGBPoint(min_gray + temp, 1.0, 0.0, 0.0);
		color->AddRGBPoint(min_gray + temp * 2, 1.0, 1.0, 0.0);
		color->AddRGBPoint(min_gray + temp * 3, 0.0, 1.0, 0.0);
		color->AddRGBPoint(min_gray + temp * 4, 0.0, 1.0, 1.0);
		color->AddRGBPoint(min_gray + temp * 5, 0.0, 0.0, 1.0);
		color->AddRGBPoint(min_gray + temp * 6, 1.0, 0.0, 1.0);
		color->AddRGBPoint(min_gray + temp * 7, 1.0, 1.0, 1.0);
		*/	
		float temp = (255 - min_gray) / 8;
		color->AddRGBPoint(min_gray, 0, 0.0, 0.0);
		color->AddRGBPoint(min_gray + temp, 0.125, 0.125, 0.125);
		color->AddRGBPoint(min_gray + temp * 2, 0.25, 0.25, 0.25);
		color->AddRGBPoint(min_gray + temp * 3, 0.375, 0.375, 0.375);
		color->AddRGBPoint(min_gray + temp * 4, 0.5, 0.5, 0.5);
		color->AddRGBPoint(min_gray + temp * 5, 0.625, 0.625, 0.0);
		color->AddRGBPoint(min_gray + temp * 6, 0.75, 0.75, 0.75);
		color->AddRGBPoint(min_gray + temp * 7, 0.875, 0.875, 0.875);
		color->AddRGBPoint(min_gray + temp * 8, 1.0, 1.0, 1.0);

		volumeProperty->SetColor(color);

		volume = vtkSmartPointer<vtkVolume>::New();
		volume->SetMapper(volumeMapper);
		volume->SetProperty(volumeProperty);

		renderer = vtkSmartPointer<vtkRenderer>::New();
		renderer->SetBackground(0.0, 0.0, 0.0);
		renderer->AddVolume(volume);
		//ren->ResetCamera();

		auto window = widget.GetRenderWindow();

		window->AddRenderer(renderer);
		window->Render();

		interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		interactor->SetRenderWindow(window);

		auto style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
		interactor->SetInteractorStyle(style);

		interactor->Initialize();

		//坐标轴
		axes = vtkSmartPointer<vtkAxesActor>::New();
		axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(1, 0, 0);
		axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0, 2, 0);
		axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0, 0, 3);
		axes->SetConeRadius(0.3);	axes->SetConeResolution(20);

		widget1 = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
		widget1->SetOutlineColor(0.9300, 0.5700, 0.1300);
		widget1->SetOrientationMarker(axes);
		widget1->SetInteractor(interactor);
		widget1->SetViewport(0.0, 0.0, 0.25, 0.25);
		widget1->SetEnabled(1);
		widget1->InteractiveOff();
	}
	else {
		return;
	}
}

void FDTDModel::OpenExcFile(void) {
	QString filename = QFileDialog::getOpenFileName(this,
		tr("Open the file"),
		"",
		tr("*.txt"));
	if (!filename.isEmpty()) {
		std::string fileStr = filename.toStdString();
		excfile = fileStr;
		FieldLJ exc;
		exc.readData(fileStr);
		exc.updateData();
		FieldLJ out;
		out.readData("./outAperture.txt");
		out.updateData();


		if (excisopen) {
			renderer->RemoveActor(fieldInActor);
			renderer->RemoveActor(fieldOutActor);
		}
		fieldInActor = exc.getActor();
		fieldOutActor = out.getActor();
		renderer->AddActor(fieldInActor);
		renderer->AddActor(fieldOutActor);

		auto window = widget.GetRenderWindow();

		//window->AddRenderer(renderer);
		window->Render();
		excisopen = true;

		if (modelisopen&excisopen)
			FDTDModel::SetFDTD->setEnabled(true);

	}
	else {
		return;
	}


}

void FDTDModel::transparencyChanged(int v) {
	transparencyLineEdit->setText(QString::number(v / 100.0f*0.1+0.9));
	transparency = v / 100.0f*0.1+0.9;
	update();
}

void FDTDModel::update(void) {
	vtkSmartPointer<vtkVolumeProperty> volumeProperty
		= vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->SetInterpolationTypeToLinear();
	volumeProperty->ShadeOn();  //打开或者关闭阴影测试
	volumeProperty->SetAmbient(1.0);
	volumeProperty->SetDiffuse(1.0);
	volumeProperty->SetSpecular(1.0);

	vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity =
		vtkSmartPointer<vtkPiecewiseFunction>::New();
	compositeOpacity->ClampingOff();
	compositeOpacity->AddPoint(min_gray, 0);
	compositeOpacity->AddPoint(max_gray, transparency);

	volumeProperty->SetScalarOpacity(compositeOpacity); //设置不透明度传输函数


	vtkSmartPointer<vtkColorTransferFunction> color =
		vtkSmartPointer<vtkColorTransferFunction>::New();
	float temp = (255 - min_gray) / 8;
	color->AddRGBPoint(min_gray, 0, 0.0, 0.0);
	color->AddRGBPoint(min_gray + temp, 0.125, 0.125, 0.125);
	color->AddRGBPoint(min_gray + temp * 2, 0.25, 0.25, 0.25);
	color->AddRGBPoint(min_gray + temp * 3, 0.375, 0.375, 0.375);
	color->AddRGBPoint(min_gray + temp * 4, 0.5, 0.5, 0.5);
	color->AddRGBPoint(min_gray + temp * 5, 0.625, 0.625, 0.0);
	color->AddRGBPoint(min_gray + temp * 6, 0.75, 0.75, 0.75);
	color->AddRGBPoint(min_gray + temp * 7, 0.875, 0.875, 0.875);
	color->AddRGBPoint(min_gray + temp * 8, 1.0, 1.0, 1.0);
	volumeProperty->SetColor(color);
	
	volume->SetMapper(volumeMapper);
	volume->SetProperty(volumeProperty);
	renderer->RemoveVolume(volume);
	renderer->AddVolume(volume);
	/*
	if (!Field3D) {
		volume->SetMapper(volumeMapper);
		volume->SetProperty(volumeProperty);
		renderer->RemoveVolume(volume);
		renderer->AddVolume(volume);

	}
	else {
		Fieldvolume->SetMapper(FieldvolumeMapper);
		Fieldvolume->SetProperty(volumeProperty);
		renderer->RemoveVolume(Fieldvolume);
		renderer->AddVolume(Fieldvolume);
	}
	*/
	

	//ren->ResetCamera();

	auto window = widget.GetRenderWindow();

	//window->AddRenderer(renderer);
	window->Render();


	interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(window);

	auto style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	interactor->SetInteractorStyle(style);
	interactor->Initialize();


	//坐标轴
	axes = vtkSmartPointer<vtkAxesActor>::New();
	axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(1, 0, 0);
	axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0, 2, 0);
	axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0, 0, 3);
	axes->SetConeRadius(0.3);	axes->SetConeResolution(20);

	widget1 = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
	widget1->SetOutlineColor(0.9300, 0.5700, 0.1300);
	widget1->SetOrientationMarker(axes);
	widget1->SetInteractor(interactor);
	widget1->SetViewport(0.0, 0.0, 0.25, 0.25);
	widget1->SetEnabled(1);
	widget1->InteractiveOff();
}

void FDTDModel::RecieveFreq(double _freq) {
	FDTDModel::frequency = _freq;
	FDTDModel::FreqEdit->setText(QString::number(frequency / 1.0e9, 10, 2));
}

void FDTDModel::SetFDTDprocess(void) {
	bool ok = false;
	frequency = FreqEdit->text().toDouble(&ok);
	frequency = frequency*1.0e9;
	OMPNumber = OMPEdit->text().toInt(&ok);
	if (!FDTD) {
		FDTD = new QTFDTDThread;
	}
	FDTD->setModelFile(modelfile);
	FDTD->setExcFile(excfile);
	FDTD->setComputation(frequency,OMPNumber,10,0,0);
	FDTD->setNewMode(true);
	CalculateFDTD->setEnabled(true);
}

void FDTDModel::RunFDTDprocess(void) {
	ButtonsOff();

	//renderer->RemoveVolume(volume);
	//auto window = widget.GetRenderWindow();
	//window->Render();
	//updateFDTDField();

	FDTD->start();
	progressBar->setHidden(false);

}

void FDTDModel::ButtonsOff(void) {

	OpenModel->setEnabled(false);
	OpenExc->setEnabled(false);
	SetFreList->setEnabled(false);
	SetFDTD->setEnabled(false);
	CalculateFDTD->setEnabled(false);
	SetAperture->setEnabled(false);
	CalAperture->setEnabled(false);
	//transparencyScrollBar->setEnabled(false);
}
void FDTDModel::ButtonsOn(void) {

	OpenModel->setEnabled(true);
	OpenExc->setEnabled(true);
	SetFreList->setEnabled(true);
	SetFDTD->setEnabled(true);
	CalculateFDTD->setEnabled(true);
	SetAperture->setEnabled(true);
	CalAperture->setEnabled(true);
	transparencyScrollBar->setEnabled(true);
}

void FDTDModel::recieveInt(int _int) {

}
void FDTDModel::recieveFloat(float _in) {
	progressBar->setValue(_in);
	//if (!fdtdcalculated) updateFDTDField();
	//可以考虑读取场并且显示！！！！
}
void FDTDModel::updateFDTDField(void) {
	renderer->RemoveVolume(volume);
	//if (Field3D) renderer->RemoveVolume(Fieldvolume);
	transparencyScrollBar->setEnabled(true);
	FILE* Fieldin;
	Fieldin = fopen("./EField3D.dat", "rb");
	int readi; float readf;
	fread(&readi, sizeof(int), 1, Fieldin);	FDTDModel::Nx = readi;
	fread(&readi, sizeof(int), 1, Fieldin);	FDTDModel::Ny = readi;
	fread(&readi, sizeof(int), 1, Fieldin);	FDTDModel::Nz = readi;

	fread(&readf, sizeof(float), 1, Fieldin); FDTDModel::dx = readf;
	fread(&readf, sizeof(float), 1, Fieldin); FDTDModel::dy = readf;
	fread(&readf, sizeof(float), 1, Fieldin); FDTDModel::dz = readf;

	fread(&readf, sizeof(float), 1, Fieldin); FDTDModel::cx = readf;
	fread(&readf, sizeof(float), 1, Fieldin); FDTDModel::cy = readf;
	fread(&readf, sizeof(float), 1, Fieldin); FDTDModel::cz = readf;

	//改变参数显示
	FDTDModel::NxEdit->setText(QString::number(Nx, 10, 0));
	FDTDModel::NyEdit->setText(QString::number(Ny, 10, 0));
	FDTDModel::NzEdit->setText(QString::number(Nz, 10, 0));
	FDTDModel::dxEdit->setText(QString::number(dx, 10, 6));
	FDTDModel::dyEdit->setText(QString::number(dy, 10, 6));
	FDTDModel::dzEdit->setText(QString::number(dz, 10, 6));
	//读取Field
	FDTDField.resize(FDTDModel::Nx*FDTDModel::Ny*FDTDModel::Nz);
	for (int k = 0; k < FDTDModel::Nz; k++) {
		for (int j = 0; j < FDTDModel::Ny; j++) {
			for (int i = 0; i <FDTDModel::Nx; i++) {

				fread(&readf, sizeof(float), 1, Fieldin);
				FDTDField[i + j*Nx + k*Nx*Ny] = readf;
			}
		}
	}
	fclose(Fieldin);

	double max = 0;
	double min = 0;
	for (int i = 0; i <Nx*Ny*Nz; i++) {
		if (FDTDField[i] > max) max = FDTDField[i];
		if (FDTDField[i] < min) min = FDTDField[i];
	}
	for (int i = 0; i <Nx*Ny*Nz; i++) {
		FDTDField[i] = (FDTDField[i]-min)/(max-min);
	}
	max = 0;
	min = 0;
	for (int i = 0; i <Nx*Ny*Nz; i++) {
		if (FDTDField[i] > max) max = FDTDField[i];
		if (FDTDField[i] < min) min = FDTDField[i];
	}


	auto img = vtkSmartPointer<vtkImageData>::New();
	auto info = vtkSmartPointer<vtkInformation>::New();
	//auto info3 = img->GetInformation();//此处也可以从img里面获取information，如
	//果information为空，则会创建一个新的，相当于调用New之后再SetInformation
	img->SetDimensions(FDTDModel::Nx, FDTDModel::Ny, FDTDModel::Nz);
	img->SetSpacing(FDTDModel::dx, FDTDModel::dy, FDTDModel::dz);
	//注意是Origin哦~
	img->SetOrigin(FDTDModel::cx - dx*Nx*0.5, FDTDModel::cy - dy*Ny*0.5, FDTDModel::cz - dz*Nz*0.5);
	img->SetScalarType(VTK_UNSIGNED_CHAR, info);
	img->SetNumberOfScalarComponents(1, info);
	img->AllocateScalars(info);

	unsigned char *ptr = (unsigned char *)img->GetScalarPointer();
	for (int i = 0; i < Nx*Ny*Nz; i++) {
		*ptr++ = float(FDTDField[i]*255) ;
	}

	volumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	volumeMapper->SetInputData(img);
	//FieldvolumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	//FieldvolumeMapper->SetInputData(img);

	vtkSmartPointer<vtkVolumeProperty> volumeProperty
		= vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->SetInterpolationTypeToLinear();
	//volumeProperty->ShadeOn();  //打开或者关闭阴影测试
	volumeProperty->SetAmbient(0.6);
	volumeProperty->SetDiffuse(0.6);
	volumeProperty->SetSpecular(0.6);

	vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity =
		vtkSmartPointer<vtkPiecewiseFunction>::New();
	compositeOpacity->ClampingOn();
	compositeOpacity->AddPoint(min_gray, 0);
	compositeOpacity->AddPoint(max_gray, transparency);

	volumeProperty->SetScalarOpacity(compositeOpacity); //设置不透明度传输函数


	vtkSmartPointer<vtkColorTransferFunction> color =
		vtkSmartPointer<vtkColorTransferFunction>::New();
	/*
	float temp = (255 - min_gray);
	color->AddRGBPoint(min_gray, 1.0, 1.0, 1.0);
	color->AddRGBPoint(min_gray + temp*0.05, 0.875, 0.875, 0.875);
	color->AddRGBPoint(min_gray + temp * 0.1, 0.75, 0.75, 0.75);
	color->AddRGBPoint(min_gray + temp * 0.2, 0.625, 0.625, 0.625);
	color->AddRGBPoint(min_gray + temp * 0.3, 0.5, 0.5, 0.5);
	color->AddRGBPoint(min_gray + temp * 0.4, 0.375, 0.375, 0.375);
	color->AddRGBPoint(min_gray + temp * 0.5, 0.25, 0.25, 0.25);
	color->AddRGBPoint(min_gray + temp * 0.7, 0.125, 0.125, 0.125);
	color->AddRGBPoint(min_gray + temp * 0.9, 0.0, 0.0, 0.0);
	*/
	float temp = (255 - min_gray);
	color->AddRGBPoint(min_gray, 0.0, 0.0, 0.0);
	color->AddRGBPoint(min_gray + temp*0.01, 1.0, 0.0, 0.0);
	color->AddRGBPoint(min_gray + temp*0.03, 1.0, 1.0, 0.0);
	color->AddRGBPoint(min_gray + temp*0.05, 0.0, 1.0, 0.0);
	color->AddRGBPoint(min_gray + temp*0.1, 0.0, 1.0, 1.0);
	color->AddRGBPoint(min_gray + temp*0.2, 0.0, 0.0, 1.0);
	color->AddRGBPoint(min_gray + temp*0.5, 1.0, 0.0, 1.0);
	color->AddRGBPoint(min_gray + temp*0.75, 1.0, 1.0, 1.0);
	
	volumeProperty->SetColor(color);

	volume = vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(volumeMapper);
	volume->SetProperty(volumeProperty);
	renderer->AddVolume(volume);

	//Fieldvolume = vtkSmartPointer<vtkVolume>::New();
	//Fieldvolume->SetMapper(FieldvolumeMapper);
	//Fieldvolume->SetProperty(volumeProperty);
	//renderer->AddVolume(Fieldvolume);
	auto window = widget.GetRenderWindow();
	window->AddRenderer(renderer);
	window->Render();

	Field3D = true;

}