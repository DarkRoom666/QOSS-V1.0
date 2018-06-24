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
#include <map>

class Mirror;
class LightShow;
class Radiator;
class Field;
class WaveguideRadiator;

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

	void setNameProject(const string & name) { nameProject = name; }

	void setUnit(double _unit) { unit = _unit; }
	double getUnit() const { return unit; }

	void setPattern(int _pattern) { pattern = _pattern; isModifiedFlag = true; }
	int getPattern() const { return pattern; }

	void setMirror(int index, Mirror* _mirror);
	Mirror* getMirror(int index) const;

	// ����Ĭ�ϵľ���
	void createDefaultMirror();
	void createModelMirror();
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

	// ����Ĭ�ϵĹ��ߣ��и�����Դ���ɣ�
	void createDefaultLigthShow();
	shared_ptr<LightShow> getDefaultLightShow() { return defaultLigthShow; }

	// ���ɷ�����
	void createRadiator();
	shared_ptr<Radiator> getRadiator() { return radiator; }

	// ���ɲ���������
	void setWaveguideRadiator(const shared_ptr<WaveguideRadiator> &);
	shared_ptr<WaveguideRadiator> getWaveguideRadiator() { return waveguideRadiator; }

	// Դ
	void setSourceField(Field*);
	Field* getSourceField() const; 

	Field* calculateByPVVA(double fre, double dis = 0.5, int N = 2);
	int addField(Field*);
	Field* getFieldByNum(int) const;

	shared_ptr<Field> getPhsCorField() const { return phsCorField; }
	void calcPhsCorField();

	bool getIsNeedCalcPhsCorFlag() const { return isNeedCalcPhsCorFlag; }

	int save(const string & dir);

	//��ע�⡿ ����open�ɹ�����޸�ԭ�е�ָ��
	static int open(const string & dir);

	void clear();

	void mesh(double ds);

private:
	static MyData * _myData;
	bool isModifiedFlag;  // ��־���������Ƿ��޸�
	bool isNeedCalcPhsCorFlag;  // ��־��λ���������Ƿ��޸�

	int numOfMirrors;
	double frequency;
	int pattern; // 0����Lower order  1����higher order 2����waveguide
	double unit;

	string nameProject;

	//����
	std::vector<Mirror*> mirrors;

	shared_ptr<LimitBox> limitBox;
	shared_ptr<calculation::SourceModeGeneration> source;
	shared_ptr<calculation::MirrorPosition> mirrorPosition;

	shared_ptr<LightShow> defaultLigthShow;
	shared_ptr<Radiator> radiator;
	shared_ptr<WaveguideRadiator> waveguideRadiator; // only in Waveguide 

	// ����������ĳ� ��0��ΪԴ��ֻ��ӵ��һ��Դ
	std::map<int, Field*> fieldMap;
	int fieldNum;

	shared_ptr<Field> phsCorField; //����������λ���������볡


};



#endif // MYDATA_H