/*
*	created by liyun 2017/11/23
*   function ����3D ������ʾ
*   version 1.0
*/

#ifndef LIGHTSHOW_H  
#define LIGHTSHOW_H

#include <vtkLine.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <list>

class LightShow
{
public:
	LightShow(int num);
	~LightShow();

	void updateData();

	vtkSmartPointer<vtkActor> getActor() const;

private:

	std::list<vtkSmartPointer<vtkActor>> actors;

	bool isTransparent; // �Ƿ���ʾ
	vtkSmartPointer<vtkProperty> property;

};

#endif //LIGHTSHOW_H