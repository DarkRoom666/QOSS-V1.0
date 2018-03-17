/*
*	created by liyun 2018/3/17
*   function ��λ������
*   version 1.0
*/

#ifndef STLMirror_H  
#define STLMirror_H

#include "Mirror.h"
#include <vector>
#include <string>

class STLMirror : public Mirror
{
public:
	STLMirror(); // ���̳���

				 // ����Ĭ�ϵ�ƽ�澵 û�д������
	STLMirror(const GraphTrans & _graphTrans);
	// ����ƽ�澵 �������
	STLMirror(const GraphTrans & _graphTrans,
		const std::vector<double> parameter);
	virtual ~STLMirror();

	virtual void calPolyData(double ds = 0);

	virtual void updateData();

	virtual QTreeWidgetItem* getTree();

	void setNameFile(const std::string & file);

	void readData();

private:
	std::string fileName;
	vtkSmartPointer<vtkPolyData> polyDataSTL;

};



#endif // STLMirror_H
