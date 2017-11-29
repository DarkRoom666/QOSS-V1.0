/*
*	created by liyun 2017/11/23
*   function ����3D ������ʾ ������ߵ����
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
#include <vector>
#include "../util/Vector3.h"
#include "../VTK/include/Mirror.h"
#include "../Calculation/SourceModeGeneration.h"
#include <memory>

class LightShow
{
public:
	LightShow(std::vector<Mirror*>, int);
	~LightShow();

	void updateData();
	void calRayActor();

	std::list<vtkSmartPointer<vtkActor>> getActors() const;

	void createStartPointBySource(shared_ptr<calculation::SourceModeGeneration>);

private:

	std::list<vtkSmartPointer<vtkActor>> actors;

	bool isTransparent; // �Ƿ���ʾ
	vtkSmartPointer<vtkProperty> property;

	int phiNum; 
	int cylinderNum;
	std::vector <std::vector <Vector3>> rayPosition;
	std::vector <std::vector <Vector3>> rayVector;
	std::vector<Mirror*> mirrors;
};

#endif //LIGHTSHOW_H