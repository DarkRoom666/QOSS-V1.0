/*
*	created by liyun 2017/10/23
*   function ���ӵ������ 
*   version 1.0
*/

#ifndef MIRROR_H  
#define MIRROR_H

#include "util/GraphTrans.h"
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>

#include <vector>
#include <../util/Vector3.h>

enum MirrorsType
{
	PLANEMIRROR = 0,
	QUADRICSURFACE
};

class actor;
class Mirror 
{
public:
	Mirror();
	virtual ~Mirror();

	MirrorsType getMirrorsType() const;

	// cal(����) and get actor 
	vtkSmartPointer<vtkActor> getActor() const;
	virtual void calActor();

	// cal(����) and get polyData 
	vtkSmartPointer<vtkPolyData> getPolyData() const;
	virtual void calPolyData(double ds = 0) = 0;

	virtual void updateData() = 0;

	// ���������߼��㽻��ͷ���
	virtual void calcReflectedRay(const Vector3&, Vector3&, Vector3&) = 0;


protected:
	GraphTrans graphTrans; // ��תƽ�Ʋ���
	vector<double> data;
	MirrorsType type;

	//����ÿ��ģ�͵���ʾ���ʷ�����
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkPolyData> polyData;

	bool isTransparent; // �Ƿ���ʾ
	vtkSmartPointer<vtkProperty> property;

private:

};



#endif // MIRROR_H