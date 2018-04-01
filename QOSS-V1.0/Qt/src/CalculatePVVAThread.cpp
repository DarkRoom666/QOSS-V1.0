#include "../include/CalculatePVVAThread.h"
#include <../Calculation/FDTDRadiator.h>
#include "../Calculation/PVVA.h"

using namespace userInterface;
using namespace calculation;
void CalculatePVVAThread::run()
{
	emit sendMainValue(0);
	MyData *myData = MyData::getInstance();
	PVVA pvva;
	pvva.SetReturnFloat(CalculatePVVAThread::setSlaverValue, this);
	// ���õ�λ
	pvva.setUnit(1);
	// ����Ƶ��
	pvva.setFre(fre);
	// ����Դ�������ڴ�
	pvva.setSource(myData->getSourceField());
	//int N = 2;
	for (int i = 1; i <= numMirror; ++i)
	{
		emit sendMainValue(i);
		pvva.setMirror(myData->getMirrorByNum(i));
		pvva.CalZ0Theta();
		pvva.Reflect();
		pvva.InterVal();		
	}
	emit sendMainValue(numMirror + 1);
	pvva.Result(dis);
	Field *tempField = new Field;
	pvva.getField(tempField);
	emit sendField(tempField);
	emit sendMainValue(numMirror + 2);
}


void CalculatePVVAThread::setSlaverValue(float val, void *user)
{
	((CalculatePVVAThread*)user)->sendSlaverValue(val);
}

void CalculatePVVAThread::killFDTD()
{
	deleteLater();
}
