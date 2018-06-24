#ifndef QTFDTDThread_H
#define QTFDTDThread_H
#include "../../FDTDRadiator/FDTDRadiator.h"
#include <QThread>
//需要调整调整了

namespace  userInterface {
	class QTFDTDThread : public QThread
	{
		Q_OBJECT

	public:
		QTFDTDThread()
		{
			FDTDradiator = new FDTDRadiator;
			newmode = false;
		}
		//Add
		void SetFDTDThread(double _frequency, double _radius, int _Ns, int _m, int _n, double _lcut, double _prodis, double _aperlen, int _Na) {
			this->frequency = _frequency;
			this->radius = _radius;
			this->Ns = _Ns;
			this->m = _m;
			this->n = _n;
			this->lcut = _lcut;
			this->prodis = _prodis;
			this->aperlen = _aperlen;
			this->Na = _Na;
			newmode = false;
		}

		~QTFDTDThread() {
			FDTDradiator->~FDTDRadiator();
		}
		void setNewMode(bool _in);
		void GetProPowerRatio(double& _PowerRatio);
		void setModelFile(string _filename);
		void setExcFile(string _filename);
		void setComputation(double _freq, int _ompNum, int _N_spa, int _timemode, int _huygensmode);
		//void setComputationMode(int _timemode, int _huygensmode);
		void setFDTDcalculated(bool _in);
		static void setMainValue(int, void*);
		static void setSlaverValue(float, void*);

	signals:
		void sendMainValue(int);
	signals:
		void sendSlaverValue(float);

		public slots:

		void killFDTD();

	protected:
		void run();

	public:
		FDTDRadiator* FDTDradiator;
	private:

		double frequency;
		double radius;		int Ns;
		//mode
		int m;
		int n;
		//Cut Height
		double lcut;
		double prodis;
		double aperlen;		int Na;
		double PowerRatio;

		std::string modelfile;
		std::string excfile;

		int OmpNum;
		int N_spa;	//10 by default
		int timemode; //0 single frequency computation; 1 multi-frequency computation
		int huygensmode; //0 five faces huygens Box; 1. cylinder faces to be established
	

		bool newmode;
		bool FDTDcalculated;


	};

}
#endif // QTFDTDThread_H