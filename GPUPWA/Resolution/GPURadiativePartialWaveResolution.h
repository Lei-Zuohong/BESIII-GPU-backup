#ifndef GPURADIATIVEPARTIALWAVERESOLUTION_H_
#define GPURADIATIVEPARTIALWAVERESOLUTION_H_

#include "GPUBasicPartialWaveResolution.h"
#include "GPUPropagatorResolution.h"
#include "GPUComputedPropagator.h"
#include "GPUGPerpStreamContainer.h"

class GPURadiativePartialWaveResolution: public GPUBasicPartialWaveResolution, public GPUGPerpStreamContainer
{
	public:
		GPURadiativePartialWaveResolution(GPUStreamTensor2 & _orbital, GPUPropagatorResolution & _prop1, char * _name);
		GPURadiativePartialWaveResolution(GPUStreamTensor2 & _orbital, GPUPropagatorResolution & _prop1, GPUPropagatorResolution & _prop2, char * _name);
		virtual ~GPURadiativePartialWaveResolution(void);
        
		unsigned int GetNProp() { return mNProp; };
		GPUPropagatorResolution & GetPropagator1() { return mProp1; };
		GPUPropagatorResolution & GetPropagator2() { return mProp2; };
		GPUStreamTensor2 & GetTensor() { return mOrbital; };
        /// Core functionality: Contract with another wave to a complex number stream
        GPUDataStream<float2> * Contract(GPUPartialWave * wave2, unsigned int index, unsigned int block);
        void IncreaseUsecount(unsigned int index) {
			mOrbital.IncreaseUsecount(index);
			GPUBasicPartialWaveResolution::IncreaseUsecount(index);
			mg_perp_stream->IncreaseUsecount(index);
		}

	protected:
		unsigned int mNProp;
		GPUPropagatorResolution & mProp1;
		GPUPropagatorResolution & mProp2;
		GPUStreamTensor2 & mOrbital;
};

#endif
