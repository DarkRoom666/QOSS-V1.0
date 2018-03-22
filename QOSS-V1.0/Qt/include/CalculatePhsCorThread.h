#ifndef CalculatePhsCorThread_H
#define CalculatePhsCorThread_H

#include "../VTK/include/Field.h"
#include <QThread>

class Mirror;
namespace  userInterface {
	class CalculatePhsCorThread : public QThread
	{
		Q_OBJECT

	public:
		CalculatePhsCorThread() { isNeedMesh = false; };
		//~CalculatePhsCorThread();

		void setDs_Length(int dsIndex, double length);

	signals:
		void sendMainValue(int);
	signals:
		void sendSlaverValue(int);
	signals:
		void sendMirror(Mirror*);
	signals:
		void error(int);

	public slots:

		void killPhsCor();

	protected:
		void run();

	private:
		Field* field;
		bool isNeedMesh;

		double ds;
		double length;
		int dsIndex; // �û�ѡ����ʷ־��� 0�� 1��׼ 2ϸ
	};

}
#endif // CalculatePhsCorThread_H