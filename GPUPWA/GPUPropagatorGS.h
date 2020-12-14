/// \file GPUPropagatorGS.h
#ifndef GPUPROPAGATORGS_H_
#define GPUPROPAGATORGS_H_

#include "GPUBasicPropagator.h"
#include "GPUPropagatorDerivative.h"
#include "ResCfg.h"

/// Implementation of a Breit-Wigner Propagator

/** Stream based Relativistic Breit-Wigner propagator. The output stream of float2s, representing complex
numbers is calculated on the first call of the () operator. Any change in the mass or width
parameters invalidates the cached stream.

\f$ BW(m_{x}^{2}) = \frac{1}{m_{r}^{2} -m_{x}^{2} + i m_{r} \Gamma_{r}} \f$
**/

class GPUPropagatorGS : public GPUBasicPropagator, public GPUPropagatorDerivative
{
public:
	/// Constructor
	GPUPropagatorGS(char *name, ///< name of the propagator
					GPUStreamScalar &_mx2,
					float _decay_m1, /// <! mass of the 1st decay particle
					float _decay_m2, /// <! mass of the 2nd decay particle
					float _r,		 /// <! effective barrier radius
					int _l			 /// <! angular momentum
	);								 ///< data stream with the squared invariant mass of the daughters
	/// Destructor
	virtual ~GPUPropagatorGS(void);
	/// Bracket operator performs the calculation or returns a cached value
	virtual GPUDataStream<float2> *operator()(int index = 0, int block = 0);
	/// Get a stream of the derivative of the propagator with regards to the mass
	virtual GPUDataStream<float2> *GetMassDerivative(int index, int block)
	{
		assert(0);
		return 0;
	};
	/// Get a stream of the derivative of the propagator with regards to the width
	virtual GPUDataStream<float2> *GetWidthDerivative(int index, int block)
	{
		assert(0);
		return 0;
	};
	/// Get a stream of the derivative of the |propagator|^2 with regards to the mass
	virtual GPUDataStream<float> *GetAbsMassDerivative(int index, int block)
	{
		assert(0);
		return 0;
	};
	/// Get a stream of the derivative of the |propagator|^2 with regards to the width
	virtual GPUDataStream<float> *GetAbsWidthDerivative(int index, int block)
	{
		assert(0);
		return 0;
	};

protected:
	/// Barrier radius
	const float m_r;
	/// Angular momentum
	const int m_l;
	/// Mass square of first decay particle
	const float m1_2;
	/// Mass square of second decay particle
	const float m2_2;
};

#endif /*GPUPROPAGATORGS_H_*/
