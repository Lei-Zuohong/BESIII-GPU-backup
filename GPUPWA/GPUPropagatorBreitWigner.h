/// \file GPUPropagatorBreitWigner.h

#pragma once
//#include "GPUBasicPropagator.h"
#include "GPUPropagatorResolution.h"
#include "GPUPropagatorDerivative.h"
#include "ResCfg.h"

/// Implementation of a Breit-Wigner Propagator

/** Stream based Breit-Wigner propagator. The output stream of float2s, representing complex
numbers is calculated on the first call of the () operator. Any change in the mass or width
parameters invalidates the cached stream.

\f$ BW(m_{x}^{2}) = \frac{1}{m_{r}^{2} -m_{x}^{2} + i m_{r} \Gamma_{r}} \f$
**/

class GPUPropagatorBreitWigner :
	public GPUPropagatorResolution, public GPUPropagatorDerivative
{
public:
	/// Constructor
	GPUPropagatorBreitWigner(char * name,  ///< name of the propagator
							GPUStreamScalar & _mx2);///< data stream with the squared invariant mass of the daughters
        GPUPropagatorBreitWigner(char * name,
                                                        GPUStreamScalar & _mx2,
                                                        float mean, float sigma);
        GPUPropagatorBreitWigner(char * name,
                                                        GPUStreamScalar & _mx2,
                                                        float mean1, float sigma1,
                                                        float mean2, float sigma2,
                                                        float frac);

	/// Destructor
	virtual ~GPUPropagatorBreitWigner(void);
	/// Bracket operator performs the calculation or returns a cached value
	virtual GPUDataStream<float2> * operator()(int index =0, int block=0);
	/// Get a stream of the derivative of the propagator with regards to the mass
	virtual GPUDataStream<float2> * GetMassDerivative(int index, int block);
	/// Get a stream of the derivative of the propagator with regards to the width
	virtual GPUDataStream<float2> * GetWidthDerivative(int index, int block);
	/// Get a stream of the derivative of the |propagator|^2 with regards to the mass
	virtual GPUDataStream<float> * GetAbsMassDerivative(int index, int block);
	/// Get a stream of the derivative of the |propagator|^2 with regards to the width
	virtual GPUDataStream<float> * GetAbsWidthDerivative(int index, int block);
        virtual GPUDataStream<float2> * GetShiftPropagator(unsigned int num, unsigned int index, unsigned int block);

protected:

};
