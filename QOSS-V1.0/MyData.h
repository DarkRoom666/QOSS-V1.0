/*
*	created by liyun 2017/10/23
*   function ������������
*   version 1.0
*/

#ifndef MYDATA_H  
#define MYDATA_H

#include <vector>
#include <memory>
#include "../Calculation/SourceModeGeneration.h"
#include "../VTK/include/LimitBox.h"
#include "MirrorPosition.h"

class Mirror;

class MyData //����
{
private:
	MyData();
	~MyData();

public:
	static MyData * getInstance();
	static const MyData * getConstInstance();

	bool getModifiedFlag() const { return isModifiedFlag; }

	void setNumOfMirrors(int num) { numOfMirrors = num; isModifiedFlag = true; }
	int getNumOfMirrors() const { return numOfMirrors; }

	void setFrequency(double fre) { frequency = fre; isModifiedFlag = true; }
	double getFrequency() const { return frequency; }

	void setUnit(double _unit) { unit = _unit; }
	double getUnit() const { return unit; }

	void setPattern(int _pattern) { pattern = _pattern; isModifiedFlag = true; }
	int getPattern() const { return pattern; }

	void setMirror(int index, Mirror* _mirror);
	Mirror* getMirror(int index) const;

	// ����Ĭ�ϵľ���
	void createDefaultMirror();
	Mirror* getMirrorByNum(int num) const;

	// ����Դ
	void setSource(const shared_ptr<calculation::SourceModeGeneration> &);
	shared_ptr<calculation::SourceModeGeneration> getSource() { return source; }

	// ���þ��ӳ�ʼλ��
	void setMirrorPosition(const shared_ptr<calculation::MirrorPosition> &);
	shared_ptr<calculation::MirrorPosition> getMirrorPosition() { return mirrorPosition; }

	// ���ú���
	void setLimitBox(const shared_ptr<LimitBox> &);
	shared_ptr<LimitBox> getLimitBox() { return limitBox; }


private:
	static MyData * _myData;
	bool isModifiedFlag;  // ��־���������Ƿ��޸�

	int numOfMirrors;
	double frequency;
	int pattern; // 0����Lower order  1����higher order 2����waveguide
	double unit;

	//����
	std::vector<Mirror*> mirrors;

	shared_ptr<LimitBox> limitBox;
	shared_ptr<calculation::SourceModeGeneration> source;
	shared_ptr<calculation::MirrorPosition> mirrorPosition;

};



#endif // MYDATA_H