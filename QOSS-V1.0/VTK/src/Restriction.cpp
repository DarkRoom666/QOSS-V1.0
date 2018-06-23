#include "..\include\Restriction.h"

#include <vtkCylinderSource.h>
#include <vtkCubeSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include "util/comUtil.h"
#include "../util/Constant_Var.h"
#include "..\Calculation\Matrix4D.h"
#include "..\Calculation\Vector3D.h"

Restriction::Restriction(RestrictionType _type)
{
	type = _type;
	property = vtkSmartPointer<vtkProperty>::New();
	property->SetOpacity(0.2);
	property->SetColor(1, 0, 0);
	//property->SetColor(1, 1, 0);
	actor = vtkSmartPointer<vtkActor>::New();
	data.resize(4);
	data[0] = 0.1;
	data[1] = 0.2;
	updateData();
}

Restriction::Restriction(const Restriction & right)
{
	property = vtkSmartPointer<vtkProperty>::New();
	property->SetOpacity(0.2);
	property->SetColor(1, 0, 0);
	//property->SetColor(1, 1, 0);
	actor = vtkSmartPointer<vtkActor>::New();
	data.resize(4);
	data[0] = right.data[0];
	data[1] = right.data[1];
	data[2] = right.data[2];
	data[3] = right.data[3];
	graphTrans = right.graphTrans;
	type = right.type;
	updateData();
}

Restriction::Restriction(const GraphTrans & _graphTrans)
{
	property = vtkSmartPointer<vtkProperty>::New();
	property->SetOpacity(0.2);
	property->SetColor(1,0,0);
	//property->SetColor(1, 1, 0);
	actor = vtkSmartPointer<vtkActor>::New();
	graphTrans = _graphTrans;
	data.resize(4);
	data[0] = 0.1;
	data[1] = 0.2;
	updateData();
}

Restriction::Restriction(const GraphTrans & _graphTrans, 
	const std::vector<double> parameter)
{
	property = vtkSmartPointer<vtkProperty>::New();
	property->SetOpacity(0.2);
	property->SetColor(1, 0, 0);
	//property->SetColor(1, 1, 0);
	actor = vtkSmartPointer<vtkActor>::New();

	graphTrans = _graphTrans;

	setData(parameter);
	updateData();
}

Restriction::~Restriction()
{
}

void Restriction::calPolyData(double ds)
{
	if (type == RES_CUBE)
	{
		float x = data[0];
		float z = data[1];

		vtkSmartPointer<vtkCubeSource> cube = vtkSmartPointer<vtkCubeSource>::New();
		cube->SetXLength(x);
		cube->SetYLength(x);
		cube->SetZLength(z);
		cube->Update();

		polyData = cube->GetOutput();
		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		// 用户自定义平移旋转 (先移动后旋转)
		transform->Translate(graphTrans.getTrans_x(),
			graphTrans.getTrans_y(), graphTrans.getTrans_z());
		transform->RotateWXYZ(graphTrans.getRotate_theta(), graphTrans.getRotate_x(),
			graphTrans.getRotate_y(), graphTrans.getRotate_z());
		transform->Translate(0, 0, z / 2);
		vtkSmartPointer<vtkTransformPolyDataFilter> TransFilter =
			vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		TransFilter->SetInputData(polyData);
		TransFilter->SetTransform(transform); //use vtkTransform (or maybe vtkLinearTransform)
		TransFilter->Update();
		polyData = TransFilter->GetOutput();

		vtkSmartPointer<vtkTriangleFilter> triangleFilter =
			vtkSmartPointer<vtkTriangleFilter>::New();
		triangleFilter->SetInputData(polyData);
		triangleFilter->Update();
		polyData = triangleFilter->GetOutput();
	}
	else if (type == RES_CYLINDER)
	{
		float r = data[0];
		float h = data[1];
		vtkSmartPointer<vtkCylinderSource> cylinder = vtkSmartPointer<vtkCylinderSource>::New();
		cylinder->SetHeight(h);
		cylinder->SetRadius(r);
		cylinder->SetCenter(0, h / 2, 0);
		cylinder->SetResolution(40);
		cylinder->Update();

		polyData = cylinder->GetOutput();
		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		// 用户自定义平移旋转 (先移动后旋转)
		transform->Translate(graphTrans.getTrans_x(),
			graphTrans.getTrans_y(), graphTrans.getTrans_z());
		transform->RotateWXYZ(graphTrans.getRotate_theta(), graphTrans.getRotate_x(),
			graphTrans.getRotate_y(), graphTrans.getRotate_z());
		transform->RotateWXYZ(90, 1, 0, 0);
		vtkSmartPointer<vtkTransformPolyDataFilter> TransFilter =
			vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		TransFilter->SetInputData(polyData);
		TransFilter->SetTransform(transform); //use vtkTransform (or maybe vtkLinearTransform)
		TransFilter->Update();
		polyData = TransFilter->GetOutput();

		vtkSmartPointer<vtkTriangleFilter> triangleFilter =
			vtkSmartPointer<vtkTriangleFilter>::New();
		triangleFilter->SetInputData(polyData);
		triangleFilter->Update();
		polyData = triangleFilter->GetOutput();
	}
	
}

vtkSmartPointer<vtkPolyData> Restriction::getPolyData() const
{
	return polyData;
}

void Restriction::updateData()
{
	calPolyData();
	calActor();
}

QTreeWidgetItem * Restriction::getTree()
{
	QTreeWidgetItem * tree1 = new QTreeWidgetItem;
	tree1->setText(0, "Restriction");
	tree1->setData(0, Qt::UserRole, QVariant(2));
	tree1->setData(1, Qt::UserRole, QVariant(-1));
	QTreeWidgetItem * treeFocus = new QTreeWidgetItem;
	treeFocus->setText(0, QString("Radius: ") + QString::number(data[0]));

	QTreeWidgetItem * treeYMax = new QTreeWidgetItem;
	treeYMax->setText(0, QString("Depth: ") + QString::number(data[1]));

	QTreeWidgetItem * tree2 = new QTreeWidgetItem;
	tree2->setText(0, "Transform");
	QTreeWidgetItem * treeTranslation = new QTreeWidgetItem;
	string tempTranslation = "Translation: ";
	tempTranslation.append(graphTrans.getTransString());
	treeTranslation->setText(0, tempTranslation.c_str());
	QTreeWidgetItem * treeRotation = new QTreeWidgetItem;
	string tempRotation = "Rotation: ";
	tempRotation.append(graphTrans.getRotateString());
	treeRotation->setText(0, tempRotation.c_str());
	tree2->addChild(treeTranslation);
	tree2->addChild(treeRotation);

	tree1->addChild(treeFocus);
	tree1->addChild(treeYMax);
	tree1->addChild(tree2);
	return tree1;
}

vtkSmartPointer<vtkActor> Restriction::getActor() const
{
	return actor;
}

void Restriction::calActor()
{
	vtkSmartPointer<vtkPolyDataMapper>mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(polyData);
	mapper->ScalarVisibilityOff();

	actor->SetMapper(mapper);
	actor->SetProperty(property);
}

void Restriction::setDataByNum(int index, double var)
{
	data[index] = var;
	updateData();
}

Json::Value Restriction::getDataJson() const
{
	Json::Value js;
	for (const auto &x : data)
	{
		js["Data"].append(x);
	}

	Json::Value jsGraphTrans;
	packGraphTransToJson(graphTrans, jsGraphTrans);
	js["graphTrans"] = jsGraphTrans;

	return js;
}

void Restriction::genRays(vector<Vector3>& star, Vector3& to, int type)
{
	if (this->type == RES_CYLINDER)
	{
		// 按周长采样
		star.clear();
		int n = 51;
		int m = 51;
		switch (type)
		{
		case 0:
		default:
			break;
		}
		double gapR = data[0] / (n - 1);
		double gapCyl = Pi / (m - 1);

		Vector3D RotateAsix(graphTrans.getRotate_x(),
			graphTrans.getRotate_y(),
			graphTrans.getRotate_z());
		Matrix4D R_rotatMatrix = Matrix4D::getRotateMatrix(
			graphTrans.getRotate_theta(), RotateAsix);
		Matrix4D R_translateMatrix = Matrix4D::getTranslateMatrix(
			graphTrans.getTrans_x(),
			graphTrans.getTrans_y(), graphTrans.getTrans_z());
		Matrix4D R_Matrix = R_translateMatrix * R_rotatMatrix;

		double x, y;
		double startR = 0;
		double startCyl = 0;
		double eachGapCyl = 0;

		Vector3 tempStar;

		//startR += gapR * 50;
		for (int i = 0; i < n; ++i)
		{
			if (startR < 0.0000001)
			{
				eachGapCyl = 2 * Pi;
			}
			else
			{
				eachGapCyl = data[0] / startR * gapCyl;
			}
			startCyl = 0;
			while (startCyl < 2 * Pi)
			{
				x = cos(startCyl) * startR;
				y = sin(startCyl) * startR;
				tempStar = R_Matrix * Vector3(x, y, 0);

				star.push_back(tempStar);
				startCyl += eachGapCyl;
			}
			startR += gapR;
		}

		to = R_rotatMatrix * Vector3(0, 0, 1);
	}
	else if (this->type == RES_CUBE)
	{
		// 按周长采样
		star.clear();
		int n = 51;
		int m = 51;
		switch (type)
		{
		case 0:
		default:
			break;
		}
		double gapN = data[0] / (n - 1);
		//double gapM = data[0] / (m - 1);

		Vector3D RotateAsix(graphTrans.getRotate_x(),
			graphTrans.getRotate_y(),
			graphTrans.getRotate_z());
		Matrix4D R_rotatMatrix = Matrix4D::getRotateMatrix(
			graphTrans.getRotate_theta(), RotateAsix);
		Matrix4D R_translateMatrix = Matrix4D::getTranslateMatrix(
			graphTrans.getTrans_x(),
			graphTrans.getTrans_y(), graphTrans.getTrans_z());
		Matrix4D R_Matrix = R_translateMatrix * R_rotatMatrix;

		double x, y;
		double startR = 0;
		double startCyl = 0;
		double eachGapCyl = 0;

		Vector3 tempStar;

		for(int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				x = -data[0] / 2 + i * gapN;
				y = -data[0] / 2 + j * gapN;
				tempStar = R_Matrix * Vector3(x, y, 0);

				star.push_back(tempStar);
			}
		to = R_rotatMatrix * Vector3(0, 0, 1);
	}
}

void Restriction::setType(RestrictionType _type)
{
	type = _type;
	updateData();
}
