#include "../include/Mirror.h"
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <vtkCaptionActor2D.h>
#include <vtkTextProperty.h>

#include <vtkSTLWriter.h>
#include "../Calculation/RayTracing.h"
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDelaunay2D.h>
#include <vtkTransform.h>
#include <vtkLine.h>
Mirror::Mirror()
	:restrictions(0)
{
	property = vtkSmartPointer<vtkProperty>::New();
	property->SetOpacity(1);
	property->SetColor(180.0 / 255.0, 180.0 / 255.0, 180.0 / 255.0);
	//property->SetColor(1, 1, 0);
	actor = vtkSmartPointer<vtkActor>::New();
	actorAxes = vtkSmartPointer<vtkAxesActor>::New();
	isTransparent = false;
	isShow = true;
}

Mirror::~Mirror()
{
	clearRestrictionAll();
}

MirrorsType Mirror::getMirrorsType() const
{
	return type;
}

vtkSmartPointer<vtkActor> Mirror::getActor() const
{
	return actor;
}

void Mirror::calActor()
{
	vtkSmartPointer<vtkPolyDataMapper>mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(polyData);
	mapper->ScalarVisibilityOff();

	actor->SetMapper(mapper);
	actor->SetProperty(property);
}

vtkSmartPointer<vtkPolyData> Mirror::getPolyData() const
{
	return polyData;
}

void Mirror::setSelected(bool flag = true)
{
	//actor->GetProperty()->SetOpacity(0.5);
	//actor->GetProperty()->SetColor(0, 1, 0);
	if (flag)
	{
		property->SetOpacity(0.2);
		property->SetColor(0, 0, 1);
	}
	else
	{
		property->SetOpacity(1);
		property->SetColor(180.0 / 255.0, 180.0 / 255.0, 180.0 / 255.0);
	}

}

void Mirror::addRestriction(Restriction *ptr)
{
	restrictions.push_back(ptr);
	updateData();
}

void Mirror::setRestriction(int num, Restriction *ptr)
{
	if (num > restrictions.size())
		return;
	delete restrictions[num];
	restrictions[num] = ptr;
	updateData();
}

void Mirror::removeRestriction(int num)
{
	if (num > restrictions.size())
		return;
	delete restrictions[num];
	auto it = restrictions.begin();
	restrictions.erase(it + num);
	updateData();
}

void Mirror::removeRestrictionAll()
{
	restrictions.clear();
}

void Mirror::clearRestrictionAll()
{
	for (int i = 0; i < restrictions.size(); ++i)
	{
		delete restrictions[i];
		restrictions[i] = nullptr;
	}
	restrictions.clear();
}

Restriction * Mirror::getRestriction(int num) const
{
	if (restrictions.empty())
		return nullptr;
	if (num > restrictions.size())
		return nullptr;
	return restrictions[num];
}

void Mirror::moveRestriction(Mirror * ptr)
{
	const vector<Restriction*>& tempRes = ptr->getRestrictionAll();
	clearRestrictionAll();
	restrictions.resize(tempRes.size());
	for (int i = 0; i < restrictions.size(); ++i)
	{
		restrictions[i] = tempRes[i];
	}
	ptr->removeRestrictionAll();
	updateData();
}

void Mirror::calcRestriction()
{
	calculation::RayTracing rayTracing(this);

	//vtkSmartPointer<vtkLine> p1 = vtkSmartPointer<vtkLine>::New();
	//vtkSmartPointer<vtkCellArray> pLineCell =
	//	vtkSmartPointer<vtkCellArray>::New();
	//int cout = 0;

	for (int i = 0; i < restrictions.size(); i++)
	{
		vtkSmartPointer<vtkPoints> points =
			vtkSmartPointer<vtkPoints>::New();
		vector<Vector3> starPoint;
		Vector3 toRay;
		Vector3 rec;
		Vector3 inter;
		bool isFlag;
		restrictions[i]->genRays(starPoint, toRay);
		for (const auto & x : starPoint)
		{
			rayTracing.calcReflect(x, toRay, rec, inter, isFlag);
			if (isFlag)
			{
				points->InsertNextPoint(inter.x, inter.y, inter.z);
				//points->InsertNextPoint(x.x, x.y, x.z);
				//points->InsertNextPoint(x.x + toRay.x, x.y + toRay.y, x.z + toRay.z);

				//p1->GetPointIds()->SetId(0, cout++);
				//p1->GetPointIds()->SetId(1, cout++);
				//pLineCell->InsertNextCell(p1);

			}
		}
		polyData = vtkSmartPointer<vtkPolyData>::New();
		polyData->SetPoints(points);

		//polyData = vtkSmartPointer<vtkPolyData>::New();
		//polyData->SetPoints(points); //获得网格模型中的几何数据：点集  
		//polyData->SetLines(pLineCell);

		vtkSmartPointer<vtkDelaunay2D> delaunay =
			vtkSmartPointer<vtkDelaunay2D>::New();
		delaunay->SetInputData(polyData);
		delaunay->Update();
		polyData = delaunay->GetOutput();
	}
}

void Mirror::switchIsTransparent()
{
	isTransparent = !isTransparent;
	if (!isShow)
		return;
	if (isTransparent) // 透明 
	{
		property->SetOpacity(0.3);
	}
	else
	{
		property->SetOpacity(1);
	}
}

void Mirror::switchIsShow()
{
	isShow = !isShow;
	if (isShow) // 显示 
	{
		if (isTransparent) // 透明 
		{
			property->SetOpacity(0.3);
		}
		else
		{
			property->SetOpacity(1);
		}
	}
	else
	{
		property->SetOpacity(0);
	}
}

vtkSmartPointer<vtkAxesActor> Mirror::getActorAxes() const
{
	return actorAxes;
}

void Mirror::calcActorAxes()
{
	actorAxes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(1, 0, 0);//修改X字体颜色为红色  
	actorAxes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0, 2, 0);//修改Y字体颜色为绿色  
	actorAxes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0, 0, 3);//修改Z字体颜色为蓝色  
	actorAxes->SetXAxisLabelText("U");
	actorAxes->SetYAxisLabelText("V");
	actorAxes->SetZAxisLabelText("N");
	actorAxes->SetConeRadius(0.1);
	actorAxes->SetConeResolution(20);
	actorAxes->SetTotalLength(0.5, 0.5, 0.5);

	vtkSmartPointer<vtkTransform> transform = 
		vtkSmartPointer<vtkTransform>::New();

	// 用户自定义平移旋转 (先移动后旋转)
	transform->Translate(graphTrans.getTrans_x(),
		graphTrans.getTrans_y(), graphTrans.getTrans_z());
	transform->RotateWXYZ(graphTrans.getRotate_theta(), 
		graphTrans.getRotate_x(), graphTrans.getRotate_y(),
		graphTrans.getRotate_z());

	actorAxes->SetUserTransform(transform);
}

void Mirror::saveSTL()
{
	vtkSmartPointer<vtkSTLWriter> writer =
		vtkSmartPointer<vtkSTLWriter>::New();
	calPolyData(0.01);
	writer->SetInputData(polyData);
	writer->SetFileName("test.stl");
	writer->Update();
}

Json::Value  Mirror::getDataJson(const string& dir, int index) const
{
	return Json::Value();
	// TODO: 在此处插入 return 语句
}

