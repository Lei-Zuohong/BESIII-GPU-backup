#ifndef GPUPROPAGATORRESOLUTION_H_
#define GPUPROPAGATORRESOLUTION_H_

#include <vector>
#include <utility>
#include "GPUDataStream.h"
#include "GPUStreamTensor.h"
#include "GPUBasicPropagator.h"

class GPUPropagatorResolution: public GPUBasicPropagator
{
	public:
		GPUPropagatorResolution(char * name, GPUStreamScalar & _mx2, std::vector<char **> parnames, unsigned int npar);
		GPUPropagatorResolution(char * name, GPUStreamScalar & _mx2, std::vector<char **> parnames, unsigned int npar, std::vector<std::pair<float, float> > & _resolution);
		GPUPropagatorResolution(char * name, GPUStreamScalar & _mx2, std::vector<char **> parnames, unsigned int npar, float _mean, float _sigma);
		GPUPropagatorResolution(char * name, GPUStreamScalar & _mx2, std::vector<char **> parnames, unsigned int npar, float _mean1, float _sigma1, float _mean2, float _sigma2, float _frac);
		//GPUPropagatorResolution(GPUStreamScalar & _mx2, int _np, float _sigma);
		//GPUPropagatorResolution(GPUStreamScalar & _mx2, float _m1, float _s1, float _m2, float _s2, float _frac);
		//GPUPropagatorResolution(GPUStreamScalar & _mx2, int _np, float _m1, float _s1, float _m2, float _s2, float _frac);
		//GPUPropagatorResolution(GPUStreamScalar & _mx2, int _npars, float * _resPars);
		virtual ~GPUPropagatorResolution(void);
		
		int GetNPoint() { return mNPoint; };
		float GetShift(int num); 
		float GetWeight(int num); 
		virtual GPUStreamScalar & GetMx2(){ return mS; };
		bool HasResolution() { return mHasResolution; };
		bool IdenticalMx2(GPUPropagatorResolution & PropRes2) { return (&mS) == (&PropRes2.GetMx2()); };
		
		virtual GPUDataStream<float2> * GetShiftPropagator(unsigned int num, unsigned int index, unsigned int block) = 0;
		void ResetResolution(unsigned int index);
		void ResetResolution(unsigned int index, unsigned int block);
		
		virtual GPUDataStream<float2> * SumResolution(unsigned int index, unsigned int block);
		virtual GPUDataStream<float2> * ContractResolution(GPUPropagatorResolution & prop2, unsigned int index, unsigned int block);
		
	protected:
		GPUStreamScalar & mS;
		std::vector<std::pair<float, float> > mResolution;
		std::vector<std::vector<GPUDataStream<float2> *> > * mShiftPropagator;
		int mNPoint;
		bool mHasResolution;
		bool mMassDependent;
};

#endif
