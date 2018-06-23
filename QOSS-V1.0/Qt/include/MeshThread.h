#ifndef MESHTHREAD_H
#define MESHTHREAD_H
#include "../MyData.h"
#include "../VTK/include/Field.h"
#include <QThread>

namespace  userInterface {
	class MeshThread : public QThread
	{
		Q_OBJECT

	public:
		MeshThread(int type)
			:type(type)
		{
		}
		//~MeshThread();
		static void setSlaverValue(float, void*);

	signals:
		void sendMainValue(int);
	signals:
		void sendSlaverValue(int);

	public slots:

		void killFDTD();

	protected:
		void run();

		int type;



	};

}
#endif // CalculatePVVAThread_H