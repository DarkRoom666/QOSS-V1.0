/*
*	created by liyun 2017/11/29
*   function �������������
*   version 1.0
*/

#ifndef RADIATOR_H  
#define RADIATOR_H

#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <memory>

#include "../../Calculation/SourceModeGeneration.h"


class Radiator
{
public:
	enum RadiatorType
	{
		LOWORDER = 0,
		HIGHT

	};

	Radiator();
	virtual ~Radiator();

	void setSource(shared_ptr<calculation::SourceModeGeneration>);

	RadiatorType getMirrorsType() const;

	virtual void calActorModel() = 0;
	vtkSmartPointer<vtkActor> getActorModel() const { return actorModel; }

	virtual void calActorRay() = 0;
	vtkSmartPointer<vtkActor> getActorRay() const { return actorRay; }

	virtual double getFirstMirrorHeight(double x);

	void getRay(std::vector <std::vector <Vector3>>&,
		std::vector <std::vector <Vector3>>&) const;

protected:

	RadiatorType type;

	//����ÿ��ģ�͵���ʾ
	vtkSmartPointer<vtkActor> actorModel;
	vtkSmartPointer<vtkActor> actorRay;

	bool isTransparent; // �Ƿ���ʾ
	vtkSmartPointer<vtkProperty> property;
	vtkSmartPointer<vtkProperty> propertyRay;

	shared_ptr<calculation::SourceModeGeneration> source;

	//����Ԩ��
	double theta;
	double radius;

	// ����
	int phiNum;
	int cylinderNum;
	std::vector <std::vector <Vector3>> rayPosition;
	std::vector <std::vector <Vector3>> rayVector;
};


#endif // RADIATOR_H

