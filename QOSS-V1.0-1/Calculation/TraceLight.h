/*
*	created by liyun 2017/10/20
*   function �����������������ʲô
*   version 1.0
*/

#ifndef TRACELIGHT_H  // ȫ��д
#define TRACELIGHT_H

namespace calculation { // ����ռ� ���������ظ�

	class TraceLight //�������ֺ�ͷ�ļ�������ͬ
	{
	public:
		TraceLight() {
			// ���б�����ʼ�� ָ��NULL
		};

		// ���м̳� һ��Ҫ��virtual
		~TraceLight() {
			// ���з������ڴ�ǵ��ͷ�
			// if(p) {
			//   delete p;
			//   p = NULL;
			// }
		};   

		// ���з����ڴ� �����д��ֵ��������� = 
		TraceLight(const TraceLight& traceLight) {
			
		}
		TraceLight operator = (const TraceLight& traceLight) {

		}

		/*
		* �򵥺������ܽ��� + ����˵�� �Լ����ز���˵��
		*/
		void fun() {};

	private:

		int* p;
	};

}

#endif // TRACELIGHT_H