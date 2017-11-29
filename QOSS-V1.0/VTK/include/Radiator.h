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

enum RadiatorType
{
	LOWORDER = 0
};


class Radiator
{
public:

	Radiator();
	virtual ~Radiator();

	void setSource(shared_ptr<calculation::SourceModeGeneration>);

	RadiatorType getMirrorsType() const;

	virtual void calActorModel() = 0;
	vtkSmartPointer<vtkActor> getActorModel() const { return actorModel; }

	virtual void calActorRay() = 0;

protected:

	RadiatorType type;

	//����ÿ��ģ�͵���ʾ
	vtkSmartPointer<vtkActor> actorModel;
	vtkSmartPointer<vtkActor> actorRay;

	bool isTransparent; // �Ƿ���ʾ
	vtkSmartPointer<vtkProperty> property;

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

