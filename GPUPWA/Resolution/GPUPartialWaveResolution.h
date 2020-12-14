#ifndef GPUPARTIALWAVERESOLUTION_H_
#define GPUPARTIALWAVERESOLUTION_H_

#include "GPUBasicPartialWaveResolution.h"
#include "GPUPropagatorResolution.h"
#include "GPUComputedPropagator.h"

class GPUPartialWaveResolution: public GPUBasicPartialWaveResolution
{
	public:
		GPUPartialWaveResolution(GPUStreamVector & _orbital, GPUPropagatorResolution & _prop1, char * _name);
		GPUPartialWaveResolution(GPUStreamVector & _orbital, GPUPropagatorResolution & _prop1, GPUPropagatorResolution & _prop2, char * _name);
		virtual ~GPUPartialWaveResolution(void);
        
		unsigned int GetNProp() { return mNProp; };
		GPUPropagatorResolution & GetPropagator1() { return mProp1; };
		GPUPropagatorResolution & GetPropagator2() { return mProp2; };
		GPUStreamVector & GetTensor() { return mOrbital; };
        /// Core functionality: Contract with another wave to a complex number stream
        GPUDataStream<float2> * Contract(GPUPartialWave * wave2, unsigned int index, unsigned int block);
        void IncreaseUsecount(unsigned int index) {
			mOrbital.IncreaseUsecount(index);
			GPUBasicPartialWaveResolution::IncreaseUsecount(index);
		}

	protected:
		unsigned int mNProp;
		GPUPropagatorResolution & mProp1;
		GPUPropagatorResolution & mProp2;
		GPUStreamVector & mOrbital;
};

#endif
