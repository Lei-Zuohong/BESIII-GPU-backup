#include "GPUFactorizedPartialWave.h"
#include "GPUUnFactorizedPartialWave2P.h"
#include "GPUUnFactorizedPartialWave.h"
#include <cstring>
//#include "GPUPropagatorResolution.h"


#ifdef USECPU
#include "Contractions_cpu.h"
using namespace Contractions_CPU;
#else
#include "Contractions.h"
using namespace Contractions_GPU;
#endif


template <typename T>
GPUUnFactorizedPartialWave2P<T>::GPUUnFactorizedPartialWave2P(
							        GPUPropagatorType<T> & _propagator,
						           	GPUScalarPropagator & _propagator1,
						                GPUPropagatorResolution & _propagator2,
							        char * _name):
							        GPUPartialWave(_propagator1.GetParset(),
								  _propagator1.GetNPars(),_name, _propagator1.GetNSets()),
mPropagator(_propagator),mPropagator1(_propagator1),mPropagator2(_propagator2)
{
	for(unsigned int i=2; i < _propagator1.GetNPars()+2; i++){
		mparnames[i] = _propagator1.GetParameterNames()[i-2];
		mparindices[i] = _propagator1.GetParameters()[i-2];
		//mlastvalues[i] = _propagator.GetLastValues()[i-2];
		for(unsigned int s=0; s <mnsets; s++)
			mlastvalues[s][i] = _propagator1.GetLastValues(s)[i-2];
		//cout <<  "Creating wave " << _propagator.GetLastValues()[i-2] << endl;
	}
}

template <typename T>
GPUUnFactorizedPartialWave2P<T>::~GPUUnFactorizedPartialWave2P(void)
{
}

template <>
GPUDataStream<float2> * GPUUnFactorizedPartialWave2P<float>::Contract(GPUPartialWave * wave2, unsigned int streamindex, unsigned int block){
	int length = mPropagator.GetLength(streamindex,block);
	GPUDataStream<float2> * aiaj = new GPUDataStream<float2>(mPropagator.GetList()->GetDeviceInterface(), length);

	GPUFactorizedPartialWave<float> * castwave2 = dynamic_cast<GPUFactorizedPartialWave<float> *>(wave2);
	if(castwave2){
		kernelcontractscalar_u_f(mPropagator.GetList()->GetDeviceInterface(),
							 mPropagator(streamindex,block),
							 castwave2->GetTensor()(streamindex,block),
							 castwave2->GetPropagator()(streamindex,block),
							 aiaj);
		return aiaj;
	}

	GPUUnFactorizedPartialWave2P<float> * ucastwave2 = dynamic_cast<GPUUnFactorizedPartialWave2P<float> *>(wave2);
	if(ucastwave2){
		kernelcontractscalar_u_u(mPropagator.GetList()->GetDeviceInterface(),
								mPropagator(streamindex,block),
								(ucastwave2->GetPropagator())(streamindex,block),
								aiaj);
		return aiaj;
	}

	cout << "Contraction for this combination of wave types not implemented, aborting!" << endl;
	exit(-1);

	// Dummy statement to avoid compiler warnings
	return aiaj;
}



template <>
GPUDataStream<float2> * GPUUnFactorizedPartialWave2P<float4>::Contract(GPUPartialWave * wave2, unsigned int streamindex, unsigned int block){
	int length = mPropagator.GetLength(streamindex,block);
	GPUDataStream<float2> * aiaj = new GPUDataStream<float2>(mPropagator.GetList()->GetDeviceInterface(), length);

/*	GPUFactorizedPartialWave<float4> * castwave2 = dynamic_cast<GPUFactorizedPartialWave<float4> *>(wave2);
	if(castwave2){
		kernelcontractmesons_v_s(mPropagator.GetList()->GetDeviceInterface(),
				                 castwave2->GetTensor()(streamindex,block),
				                 mPropagator.realpart(streamindex,block),
				                 mPropagator.imagpart(streamindex,block),
				                 castwave2->GetPropagator()(streamindex,block),
				                 aiaj);
		return aiaj;
	}
*/

	GPUUnFactorizedPartialWave<float4> * ucastwave2kstar = dynamic_cast<GPUUnFactorizedPartialWave<float4> *>(wave2);
	if(ucastwave2kstar){
                                GPUPropagatorResolution & prop1a = GetPropagator2();
                                GPUDataStream<float2> * sum1 = prop1a.SumResolution(streamindex, block);
				kernelcontractmesons_v_vs(mPropagator.GetList()->GetDeviceInterface(),
								mPropagator.realpart(streamindex,block),
								mPropagator.imagpart(streamindex,block),
								ucastwave2kstar->GetPropagator().realpart(streamindex,block),
							        ucastwave2kstar->GetPropagator().imagpart(streamindex,block),
                                                                sum1,
								aiaj);
                                delete sum1;

		return aiaj;
	}




	GPUUnFactorizedPartialWave2P<float4> * ucastwave2 = dynamic_cast<GPUUnFactorizedPartialWave2P<float4> *>(wave2);
	if(ucastwave2){
//                       cout<<"1111111111"<<endl;
                        GPUPropagatorResolution & prop1a = GetPropagator2();
                        GPUPropagatorResolution & prop2a = ucastwave2->GetPropagator2();
                        if (prop1a.IdenticalMx2(prop2a)) {
//                        cout<<"22222222"<<endl;
                              GPUDataStream<float2> * contract = prop1a.ContractResolution(prop2a, streamindex, block); 
					kernelcontractmesons_v_vs(mPropagator.GetList()->GetDeviceInterface(),
								mPropagator.realpart(streamindex,block),
								mPropagator.imagpart(streamindex,block),
								ucastwave2->GetPropagator().realpart(streamindex,block),
								ucastwave2->GetPropagator().imagpart(streamindex,block),
                                                                contract,
								aiaj);
                                delete contract;
                        }
/*else{
                              cout<<"333333333"<<endl; 
                                GPUDataStream<float2> * sum1 = prop1a.SumResolution(streamindex, block);
                                kernelcontractmesons_v_vs(mPropagator.GetList()->GetDeviceInterface(),
                                                                mPropagator.realpart(streamindex,block),
                                                                mPropagator.imagpart(streamindex,block),
                                                                ucastwave2->GetPropagator().realpart(streamindex,block),
                                                                ucastwave2->GetPropagator().imagpart(streamindex,block),
                                                                sum1,aiaj);
                                delete sum1;
                        }
                  */


		return aiaj;
	}

	cout << "Contraction for this combination of wave types not implemented, aborting!" << endl;
	exit(-1);

	// Dummy statement to avoid compiler warnings
	return aiaj;
}

// dummy function, overriden by derived classes
template <>
GPUDataStream<float2> * GPUUnFactorizedPartialWave2P<float44>::Contract(GPUPartialWave * wave2, unsigned int streamindex, unsigned int block){
	assert(0);
	return 0;
}

template class GPUUnFactorizedPartialWave2P<float>;
template class GPUUnFactorizedPartialWave2P<float4>;
template class GPUUnFactorizedPartialWave2P<float44>;
