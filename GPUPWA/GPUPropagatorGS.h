#ifndef GPUPROPAGATORGS_H_
#define GPUPROPAGATORGS_H_
#include "GPUBasicPropagator.h"
#include "GPUPropagatorDerivative.h"
#include "ResCfg.h"

class GPUPropagatorGS : public GPUBasicPropagator, public GPUPropagatorDerivative
{
public:
	/// Constructor
	GPUPropagatorGS(char *name,			   ///< name of the propagator
					GPUStreamScalar &_mx2, ///
					float _decay_m1,	   /// <! mass of the 1st decay particle
					float _decay_m2,	   /// <! mass of the 2nd decay particle
					float _r,			   /// <! effective barrier radius
					int _l				   /// <! angular momentum
	);									   ///< data stream with the squared invariant mass of the daughters
	virtual ~GPUPropagatorGS(void);
	virtual GPUDataStream<float2> *operator()(int index = 0, int block = 0);
	virtual GPUDataStream<float2> *GetMassDerivative(int index, int block)
	{
		assert(0);
		return 0;
	};
	virtual GPUDataStream<float2> *GetWidthDerivative(int index, int block)
	{
		assert(0);
		return 0;
	};
	virtual GPUDataStream<float> *GetAbsMassDerivative(int index, int block)
	{
		assert(0);
		return 0;
	};
	virtual GPUDataStream<float> *GetAbsWidthDerivative(int index, int block)
	{
		assert(0);
		return 0;
	};

protected:
	const float m_r;
	const int m_l;
	const float m1_2;
	const float m2_2;
};

#endif
