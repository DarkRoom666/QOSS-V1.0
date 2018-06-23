#include "../include/MeshThread.h"
#include "../MyData.h"
#include "../VTK/include/Mirror.h"
#include "../util/Constant_Var.h"
using namespace userInterface;
using namespace calculation;
void MeshThread::run()
{
	int numMirror = MyData::getInstance()->getNumOfMirrors();
	double fre = MyData::getInstance()->getFrequency();
	double lamda = C_Speed / fre;
	double ds = lamda / 6;
	switch (type)
	{
	case 0:
		ds = lamda / 3;
		break;
	case 2:
		ds = lamda / 8;
		break;
	case 1:
	default:
		break;
	}
	emit sendMainValue(0);
	
	for (int i = 0; i < numMirror; ++i)
	{
		MyData::getInstance()->getMirrorByNum(i)->genMesh(ds);
		emit sendMainValue(i);	
	}
	emit sendMainValue(numMirror + 1);
}


void MeshThread::setSlaverValue(float val, void *user)
{
	((MeshThread*)user)->sendSlaverValue(val);
}

void MeshThread::killFDTD()
{
	deleteLater();
}
