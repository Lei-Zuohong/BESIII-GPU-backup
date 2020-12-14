#include "GPUBasicPartialWaveResolution.h"
#include <cstring>
#include <vector>

GPUBasicPartialWaveResolution::GPUBasicPartialWaveResolution(char *_name, unsigned int nsets, GPUScalarPropagator &_propagator):
    GPUPartialWave(_propagator.GetParset(), _propagator.GetNPars(), _name, nsets),
	mPropagator(_propagator) {

	for(unsigned int i=2; i < _propagator.GetNPars()+2; ++i){
		mparnames[i] = _propagator.GetParameterNames()[i-2];
		mparindices[i] = _propagator.GetParameters()[i-2];
		for(unsigned int s=0; s <mnsets; ++s) {
			mlastvalues[s][i] = _propagator.GetLastValues(s)[i-2];
		}
	}
}

GPUBasicPartialWaveResolution::~GPUBasicPartialWaveResolution() {
}
