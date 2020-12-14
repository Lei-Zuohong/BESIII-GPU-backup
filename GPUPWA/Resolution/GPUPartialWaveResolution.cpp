#include "GPUPartialWaveResolution.h"

#ifdef USECPU
#include "Resolution_cpu.h"
using namespace Resolution_CPU;
#else
#include "Resolution.h"
using namespace Resolution_GPU;
#endif
		
GPUPartialWaveResolution::GPUPartialWaveResolution(GPUStreamVector & _orbital, GPUPropagatorResolution & _prop1, char * _name):
	GPUBasicPartialWaveResolution(_name, _orbital.GetNSets(), _prop1),
	mNProp(1), mProp1(_prop1), mProp2(_prop1), mOrbital(_orbital)
{
}

GPUPartialWaveResolution::GPUPartialWaveResolution(GPUStreamVector & _orbital, GPUPropagatorResolution & _prop1, GPUPropagatorResolution & _prop2, char * _name):
	GPUBasicPartialWaveResolution(_name, _orbital.GetNSets(), _prop1*_prop2),
	mNProp(2), mProp1(_prop1), mProp2(_prop2), mOrbital(_orbital)
{
	assert(!_prop1.IdenticalMx2(_prop2));
}

GPUPartialWaveResolution::~GPUPartialWaveResolution(void)
{
}
        
GPUDataStream<float2> * GPUPartialWaveResolution::Contract(GPUPartialWave * wave2, unsigned int index, unsigned int block)
{
	int length = mOrbital.GetLength(index, block);
	GPUDataStream<float2> * aiaj = new GPUDataStream<float2>(mOrbital.GetList()->GetDeviceInterface(), length);
	
	GPUPartialWaveResolution * castwave2 = dynamic_cast<GPUPartialWaveResolution *>(wave2);
	if (castwave2) {
		unsigned int nprop1 = GetNProp();
		unsigned int nprop2 = castwave2->GetNProp();

		if (nprop1 == 1 && nprop2 == 1) {
			GPUPropagatorResolution & prop1a = GetPropagator1();
			GPUPropagatorResolution & prop2a = castwave2->GetPropagator1();
			if (prop1a.IdenticalMx2(prop2a)) {
				GPUDataStream<float2> * contract = prop1a.ContractResolution(prop2a, index, block);
				kernelcontract_reswave_10(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), contract, 
						castwave2->GetTensor()(index, block),
						aiaj);
				delete contract;
			}
			else {
				GPUDataStream<float2> * sum1 = prop1a.SumResolution(index, block);
				GPUDataStream<float2> * sum2 = prop2a.SumResolution(index, block);
				kernelcontract_reswave_11(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), sum1, 
						castwave2->GetTensor()(index, block), sum2,
						aiaj);
				delete sum1;
				delete sum2;
			}
		}
		else if (nprop1 == 1 && nprop2 == 2) {
			GPUPropagatorResolution & prop1a = GetPropagator1();
			GPUPropagatorResolution & prop2a = castwave2->GetPropagator1();
			GPUPropagatorResolution & prop2b = castwave2->GetPropagator2();
			if (prop1a.IdenticalMx2(prop2a)) {
				GPUDataStream<float2> * contract = prop1a.ContractResolution(prop2a, index, block);
				GPUDataStream<float2> * sum = prop2b.SumResolution(index, block);
				kernelcontract_reswave_11(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), contract, 
						castwave2->GetTensor()(index, block), sum, 
						aiaj);
				delete contract;
				delete sum;
			}
			else if (prop1a.IdenticalMx2(prop2b)) {
				GPUDataStream<float2> * contract = prop1a.ContractResolution(prop2b, index, block);
				GPUDataStream<float2> * sum = prop2a.SumResolution(index, block);
				kernelcontract_reswave_11(mOrbital.GetList()->GetDeviceInterface(),
					mOrbital(index, block), contract, 
					castwave2->GetTensor()(index, block), sum,
					aiaj);
				delete contract;
				delete sum;
			}
			else {
				GPUDataStream<float2> * sum1a = prop1a.SumResolution(index, block);
				GPUDataStream<float2> * sum2a = prop2a.SumResolution(index, block);
				GPUDataStream<float2> * sum2b = prop2b.SumResolution(index, block);
				kernelcontract_reswave_12(mOrbital.GetList()->GetDeviceInterface(),
					mOrbital(index, block), sum1a, 
					castwave2->GetTensor()(index, block), sum2a, sum2b,
					aiaj);
				delete sum1a;
				delete sum2a;
				delete sum2b;
			}
		}
		else if (nprop1 == 2 && nprop1 == 1) {
			GPUPropagatorResolution & prop1a = GetPropagator1();
			GPUPropagatorResolution & prop1b = GetPropagator2();
			GPUPropagatorResolution & prop2a = castwave2->GetPropagator1();
			if (prop1a.IdenticalMx2(prop2a)) {
				GPUDataStream<float2> * contract = prop1a.ContractResolution(prop2a, index, block);
				GPUDataStream<float2> * sum = prop1b.SumResolution(index, block);
				kernelcontract_reswave_20(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), contract, sum, 
						castwave2->GetTensor()(index, block),
						aiaj);
				delete contract;
				delete sum;
			}
			else if (prop1b.IdenticalMx2(prop2a)) {
				GPUDataStream<float2> * contract = prop1b.ContractResolution(prop2a, index, block);
				GPUDataStream<float2> * sum = prop1a.SumResolution(index, block);
				kernelcontract_reswave_20(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), contract, sum, 
						castwave2->GetTensor()(index, block),
						aiaj);
				delete contract;
				delete sum;
			}
			else {
				GPUDataStream<float2> * sum1a = prop1a.SumResolution(index, block);
				GPUDataStream<float2> * sum1b = prop1b.SumResolution(index, block);
				GPUDataStream<float2> * sum2a = prop2a.SumResolution(index, block);
				kernelcontract_reswave_21(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), sum1a, sum1b,
						castwave2->GetTensor()(index, block), sum2a,
						aiaj);
				delete sum1a;
				delete sum1b;
				delete sum2a;
			}
		}
		else if (nprop1 == 2 && nprop2 == 2) {
			GPUPropagatorResolution & prop1a = GetPropagator1();
			GPUPropagatorResolution & prop1b = GetPropagator2();
			GPUPropagatorResolution & prop2a = castwave2->GetPropagator1();
			GPUPropagatorResolution & prop2b = castwave2->GetPropagator2();
			if (prop1a.IdenticalMx2(prop2a) && prop1b.IdenticalMx2(prop2b)) {
				GPUDataStream<float2> * contract1 = prop1a.ContractResolution(prop2a, index, block);
				GPUDataStream<float2> * contract2 = prop1b.ContractResolution(prop2b, index, block);
				kernelcontract_reswave_20(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), contract1, contract2, 
						castwave2->GetTensor()(index, block), 
						aiaj);
				delete contract1;
				delete contract2;
			}
			else if (prop1a.IdenticalMx2(prop2a) && !prop1b.IdenticalMx2(prop2b)) {
				GPUDataStream<float2> * contract = prop1a.ContractResolution(prop2a, index, block);
				GPUDataStream<float2> * sum1b = prop1b.SumResolution(index, block);
				GPUDataStream<float2> * sum2b = prop2b.SumResolution(index, block);
				kernelcontract_reswave_21(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), contract, sum1b, 
						castwave2->GetTensor()(index, block), sum2b, 
						aiaj);
				delete contract;
				delete sum1b;
				delete sum2b;
			}
			else if (prop1a.IdenticalMx2(prop2b) && prop1b.IdenticalMx2(prop2a)) {
				GPUDataStream<float2> * contract1 = prop1a.ContractResolution(prop2b, index, block);
				GPUDataStream<float2> * contract2 = prop1b.ContractResolution(prop2a, index, block);
				kernelcontract_reswave_20(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), contract1, contract2, 
						castwave2->GetTensor()(index, block), 
						aiaj);
				delete contract1;
				delete contract2;
			}
			else if (prop1a.IdenticalMx2(prop2b) && !prop1b.IdenticalMx2(prop2a)) {
				GPUDataStream<float2> * contract = prop1a.ContractResolution(prop2b, index, block);
				GPUDataStream<float2> * sum1b = prop1b.SumResolution(index, block);
				GPUDataStream<float2> * sum2a = prop2a.SumResolution(index, block);
				kernelcontract_reswave_21(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), contract, sum1b, 
						castwave2->GetTensor()(index, block), sum2a, 
						aiaj);
				delete contract;
				delete sum1b;
				delete sum2a;
			}
			else if (prop1b.IdenticalMx2(prop2a)) {
				GPUDataStream<float2> * contract = prop1b.ContractResolution(prop2a, index, block);
				GPUDataStream<float2> * sum1a = prop1a.SumResolution(index, block);
				GPUDataStream<float2> * sum2b = prop2b.SumResolution(index, block);
				kernelcontract_reswave_21(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), contract, sum1a, 
						castwave2->GetTensor()(index, block), sum2b, 
						aiaj);
				delete contract;
				delete sum1a;
				delete sum2b;
			}
			else if (prop1b.IdenticalMx2(prop2b)) {
				GPUDataStream<float2> * contract = prop1b.ContractResolution(prop2b, index, block);
				GPUDataStream<float2> * sum1a = prop1a.SumResolution(index, block);
				GPUDataStream<float2> * sum2a = prop2a.SumResolution(index, block);
				kernelcontract_reswave_21(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), contract, sum1a, 
						castwave2->GetTensor()(index, block), sum2a, 
						aiaj);
				delete contract;
				delete sum1a;
				delete sum2a;
			}
			else {
				GPUDataStream<float2> * sum1a = prop1a.SumResolution(index, block);
				GPUDataStream<float2> * sum1b = prop1b.SumResolution(index, block);
				GPUDataStream<float2> * sum2a = prop2a.SumResolution(index, block);
				GPUDataStream<float2> * sum2b = prop2b.SumResolution(index, block);
				kernelcontract_reswave_22(mOrbital.GetList()->GetDeviceInterface(),
						mOrbital(index, block), sum1a, sum1b, 
						castwave2->GetTensor()(index, block), sum2a, sum2b, 
						aiaj);
				delete sum1a;
				delete sum1b;
				delete sum2a;
				delete sum2b;
			}
		}
		else {
			std::cout << "Contraction between " << nprop1 << " propagator(s) and " << nprop2 << " propagator(s) is not supported." << std::endl;
			exit(-1);
			return aiaj;
		}
		return aiaj;
	}
	std::cout << "Contraction for this combination of wave types not implemented, aborting!" << std::endl;
	exit(-1);
	return aiaj;
}
