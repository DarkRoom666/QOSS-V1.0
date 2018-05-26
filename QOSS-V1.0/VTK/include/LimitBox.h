/*
*	created by liyun 2017/11/23
*   function ���Ƶĺ��Ӵ�С 
*   version 1.0
*/

#ifndef LIMITBOX_H
#define LIMITBOX_H

#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkjsoncpp/json/json.h>

class LimitBox
{
public:
	LimitBox(double length, double width, double height);
	~LimitBox();

	void setLength(double _length);
	double getLength() const { return length; }

	void setWidth(double _width);
	double getWidth() const { return width; }

	void setHeight(double _height);
	double getHeight() const { return height; }

	void updateData();

	vtkSmartPointer<vtkActor> getActor() const;

	double getMaxSize() const;

	void setIsTransparent(bool);

	Json::Value getDataJson() const;

	int setDataJson(const Json::Value &);

private:
	vtkSmartPointer<vtkActor> actor;

	bool isTransparent; // �Ƿ���ʾ
	vtkSmartPointer<vtkProperty> property;

	double length;
	double width;
	double height;
};

#endif // LIMITBOX_H
