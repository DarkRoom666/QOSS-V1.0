/*
*	created by liyun 2017/10/17
*   function ����QLabel ���Ʋ���2d��ʾ
*   version 1.0
*   ����ʱ������setNum() ����������setData()
*/
/*
	Ӧ�����S�ṩ�Ļ�ͼ����г���ͼ- ֻʹ�ã����޸���ã�
*/

#ifndef QPAINTFIELD_H
#define QPAINTFIELD_H

#include <QLabel>
#include <QPainter>
#include <QPen>
#include <QPoint>
#include <vector>
#include <complex>

using namespace std;

class QPaintField : public QLabel
{

	Q_OBJECT

public:

	explicit QPaintField(QWidget *parent = 0);

	void paintEvent(QPaintEvent *event);

	// ����n*m�ľ���
	void setNum(int, int);

	// ������ʾ���� ����һ��n*m�Ķ�ά����
	void setData(vector <vector <double>> &);

	// ������������ ����ƽ���� ��������ܲ����ã����Ժ��ԣ�
	void calcData(vector <vector <complex <double>>> &,
		vector <vector <complex <double>>> &);

private:

	int NumOfLength;
	int NumOfWidth;

	vector <vector <double>> data;
};

#endif //QPAINTFIELD_H
