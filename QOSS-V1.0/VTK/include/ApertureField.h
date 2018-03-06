/*
*	created by liyun 2017/12/28
*   function ���ɸ�˹���泡
*	width --> data[0]  depth --> data[1]
*	scale --> data[2]  phase --> data[3]
*	line  --> data[4]  Column--> data[5]
*   isAmPhs  --> data[6] 0Ϊ�� 1Ϊ��
*   version 1.0
*/

#ifndef APERTUREFIELD_H
#define APERTUREFIELD_H

#include <string>
#include <complex>
#include <vector>
#include <QString>

#include "Field.h"


using namespace std;

class ApertureField : public Field
{
public:
	ApertureField(const GraphTrans & _graphTrans,
		const std::vector<double> parameter);

	ApertureField(const ApertureField & _ApertureField);

	~ApertureField();

	void setFileAddress(const string & file);

	void setUnit(double _unit);

	void updateData();

	void readData();

	virtual QTreeWidgetItem* getTree();

private:

	double factor;

	string fileAddress;

	bool isAmPhs; // �Ƿ�����λ����ʵ���鲿

};

#endif
