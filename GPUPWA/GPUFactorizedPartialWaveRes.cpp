#include "GPUFactorizedPartialWaveRes.h"
#include "GPUUnFactorizedPartialWave.h"
#include <cstring>

#ifdef USECPU
#include "Contractions_cpu.h"
using namespace Contractions_CPU;
#else
#include "Contractions.h"
using namespace Contractions_GPU;
#endif

template <typename T>
GPUFactorizedPartialWaveRes<T>::GPUFactorizedPartialWaveRes(GPUStreamTensor<T> & _orbitalTensor,
								  GPUPropagatorResolution & _propagator,
								  char * _name):
								  GPUPartialWave(_propagator.GetParset(),
								  _propagator.GetNPars(),_name, _orbitalTensor.GetNSets()),
mOrbitalTensor(_orbitalTensor),
mPropagator(_propagator)
{
	for(unsigned int i=2; i < _propagator.GetNPars()+2; i++){
		mparnames[i] = _propagator.GetParameterNames()[i-2];
		mparindices[i] = _propagator.GetParameters()[i-2];
		for(unsigned int s=0; s <mnsets; s++)
			mlastvalues[s][i] = _propagator.GetLastValues(s)[i-2];
		//cout <<  "Creating wave " << _propagator.GetLastValues()[i-2] << endl;
	}
}

template <typename T>
GPUFactorizedPartialWaveRes<T>::~GPUFactorizedPartialWaveRes(void)
{
	delete[] mName;
}


template <>
GPUDataStream<float2> * GPUFactorizedPartialWaveRes<float>::Contract(GPUPartialWave * wave2, unsigned int streamindex, unsigned int block){
        assert(0);
        return 0;
}

template <>
GPUDataStream<float2> * GPUFactorizedPartialWaveRes<float4>::Contract(GPUPartialWave * wave2, unsigned int streamindex, unsigned int block){
	int length = mOrbitalTensor.GetLength(streamindex,block);
	GPUDataStream<float2> * aiaj = new GPUDataStream<float2>(mOrbitalTensor.GetList()->GetDeviceInterface(), length);

	GPUFactorizedPartialWaveRes<float4> * castwave2 = dynamic_cast<GPUFactorizedPartialWaveRes<float4> *>(wave2);
	if(castwave2){
		GPUPropagatorResolution & prop1a = GetPropagator();
		GPUPropagatorResolution & prop2a = castwave2->GetPropagator();

		if (prop1a.IdenticalMx2(prop2a)) {

			GPUDataStream<float2> * contract = prop1a.ContractResolution(prop2a, streamindex, block);

			kernelcontract_reswave_10(mOrbitalTensor.GetList()->GetDeviceInterface(),
						  mOrbitalTensor(streamindex, block), contract, 
                                                  castwave2->GetTensor()(streamindex, block),
                                                  aiaj);

			delete contract;
		}

//			ValidateCache();
//			castwave2->ValidateCache();
			return aiaj;
		
	}	

	GPUUnFactorizedPartialWave<float4> * ucastwave2 = dynamic_cast<GPUUnFactorizedPartialWave<float4> *>(wave2);
        	if(ucastwave2){

		GPUPropagatorResolution & prop1a = GetPropagator();
		GPUDataStream<float2> * sum1 = prop1a.SumResolution(streamindex, block);

		kernelcontractmesons_s_v(mOrbitalTensor.GetList()->GetDeviceInterface(),
									mOrbitalTensor(streamindex,block),
									sum1,
									(ucastwave2->GetPropagator().realpart(streamindex,block)),
									(ucastwave2->GetPropagator().imagpart(streamindex,block)),
									aiaj);
		ValidateCache();
		ucastwave2->ValidateCache();
		return aiaj;
	}

	cout << "Contraction for this combination of wave types not implemented, aborting!" << endl;
	exit(-1);

	// Dummy statement to avoid compiler warnings
	return aiaj;
}
template <>
GPUDataStream<float2> * GPUFactorizedPartialWaveRes<float44>::Contract(GPUPartialWave * wave2, unsigned int streamindex, unsigned int block){
        assert(0);
        return 0;
}


/// Dummy fuctions, overridden by derived classes but neede in linking

template class GPUFactorizedPartialWaveRes<float>;
template class GPUFactorizedPartialWaveRes<float4>;
template class GPUFactorizedPartialWaveRes<float44>;
//template class GPUPartialWave<float44>;
