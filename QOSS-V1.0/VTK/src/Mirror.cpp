#include "../include/Mirror.h"
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <vtkCaptionActor2D.h>
#include <vtkTextProperty.h>

#include <vtkSTLWriter.h>

Mirror::Mirror()
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

void Mirror::switchIsTransparent()
{
	isTransparent = !isTransparent;
	if (!isShow)
		return;
	if (isTransparent) // ͸�� 
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
	if (isShow) // ��ʾ 
	{
		if (isTransparent) // ͸�� 
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
	actorAxes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(1, 0, 0);//�޸�X������ɫΪ��ɫ  
	actorAxes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0, 2, 0);//�޸�Y������ɫΪ��ɫ  
	actorAxes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0, 0, 3);//�޸�Z������ɫΪ��ɫ  
	actorAxes->SetXAxisLabelText("U");
	actorAxes->SetYAxisLabelText("V");
	actorAxes->SetZAxisLabelText("N");
	actorAxes->SetConeRadius(0.1);
	actorAxes->SetConeResolution(20);
	actorAxes->SetTotalLength(0.5, 0.5, 0.5);

	vtkSmartPointer<vtkTransform> transform = 
		vtkSmartPointer<vtkTransform>::New();

	// �û��Զ���ƽ����ת (���ƶ�����ת)
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
	// TODO: �ڴ˴����� return ���
}

