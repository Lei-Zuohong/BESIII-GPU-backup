#include "GPUPropagatorResolution.h"

#ifdef USECPU
#include "Propagators_cpu.h"
using namespace Propagators_CPU;
#else
#include "Propagators.h"
using namespace Propagators_GPU;
#endif

#include <iostream>

// no resolution
GPUPropagatorResolution::GPUPropagatorResolution(char * name, GPUStreamScalar & _mx2, std::vector<char **> parnames, unsigned int npar):
	GPUBasicPropagator(name, _mx2, std::vector<char**>(npar,(char **)0), npar),
	mS(_mx2), mNPoint(1), mHasResolution(false), mMassDependent(false)
{
	mResolution.push_back(std::make_pair(0, 1));

	mShiftPropagator = new std::vector<std::vector<GPUDataStream<float2> *> >[mNPoint];
	for (int i = 0; i != mNPoint; ++i) {
		mShiftPropagator[i] = std::vector<std::vector<GPUDataStream<float2> *> >(mS.GetList()->GetNSets(), std::vector<GPUDataStream<float2> *>());
		for	(int j = 0; j != mS.GetList()->GetNSets(); ++j){
			for (int k = 0; k != mS.GetList()->GetNBlocks(j); ++k){
				mShiftPropagator[i][j].push_back(0);
			}
		}
	}
}

// set resolution by shift-weight pair
GPUPropagatorResolution::GPUPropagatorResolution(char * name, GPUStreamScalar & _mx2, std::vector<char **> parnames, unsigned int npar, std::vector<std::pair<float, float> > & _resolution):
	GPUBasicPropagator(name, _mx2, std::vector<char**>(npar,(char **)0), npar),
	mS(_mx2), mResolution(_resolution), mNPoint(_resolution.size()), mHasResolution(true), mMassDependent(false)
{
	assert(_resolution.size() >= 1);
	
	mShiftPropagator = new std::vector<std::vector<GPUDataStream<float2> *> >[mNPoint];
	for (int i = 0; i != mNPoint; ++i) {
		mShiftPropagator[i] = std::vector<std::vector<GPUDataStream<float2> *> >(mS.GetList()->GetNSets(), std::vector<GPUDataStream<float2> *>());
		for	(int j = 0; j != mS.GetList()->GetNSets(); ++j){
			for (int k = 0; k != mS.GetList()->GetNBlocks(j); ++k){
				mShiftPropagator[i][j].push_back(0);
			}
		}
	}
}

// gaussian resolution, 11-point gauss-hermite integration
GPUPropagatorResolution::GPUPropagatorResolution(char * name, GPUStreamScalar & _mx2, std::vector<char **> parnames, unsigned int npar, float _mean, float _sigma): 
	GPUBasicPropagator(name, _mx2, std::vector<char**>(npar,(char **)0), npar),
	mS(_mx2), mNPoint(11), mHasResolution(true), mMassDependent(false)
{
	mShiftPropagator = new std::vector<std::vector<GPUDataStream<float2> *> >[mNPoint];
	for (int i = 0; i != mNPoint; ++i) {
		mShiftPropagator[i] = std::vector<std::vector<GPUDataStream<float2> *> >(mS.GetList()->GetNSets(), std::vector<GPUDataStream<float2> *>());
		for	(int j = 0; j != mS.GetList()->GetNSets(); ++j){
			for (int k = 0; k != mS.GetList()->GetNBlocks(j); ++k){
				mShiftPropagator[i][j].push_back(0);
			}
		}
	}

	mResolution.push_back(std::make_pair(-_mean-5.188000e+00f*_sigma, 8.121848e-07f));
	mResolution.push_back(std::make_pair(-_mean-3.936166e+00f*_sigma, 1.956717e-04f));
	mResolution.push_back(std::make_pair(-_mean-2.865126e+00f*_sigma, 6.720288e-03f));
	mResolution.push_back(std::make_pair(-_mean-1.876039e+00f*_sigma, 6.613882e-02f));
	mResolution.push_back(std::make_pair(-_mean-9.288696e-01f*_sigma, 2.422404e-01f));
	mResolution.push_back(std::make_pair(-_mean+0.000000e+00f*_sigma, 3.694082e-01f));
	mResolution.push_back(std::make_pair(-_mean+9.288696e-01f*_sigma, 2.422404e-01f));
	mResolution.push_back(std::make_pair(-_mean+1.876039e+00f*_sigma, 6.613882e-02f));
	mResolution.push_back(std::make_pair(-_mean+2.865126e+00f*_sigma, 6.720288e-03f));
	mResolution.push_back(std::make_pair(-_mean+3.936166e+00f*_sigma, 1.956717e-04f));
	mResolution.push_back(std::make_pair(-_mean+5.188000e+00f*_sigma, 8.121848e-07f));
}

//// gaussian resolution
//GPUPropagatorResolution::GPUPropagatorResolution(GPUStreamScalar & _mx2, int _np, float _sigma): 
//	mS(_mx2), mNPoint(_np), mHasResolution(true), mMassDependent(false)
//{
//	assert(_np > 1);
//}

// double gaussian resolution
GPUPropagatorResolution::GPUPropagatorResolution(char * name, GPUStreamScalar & _mx2, std::vector<char **> parnames, unsigned int npar, float _mean1, float _sigma1, float _mean2, float _sigma2, float _frac):
	GPUBasicPropagator(name, _mx2, std::vector<char**>(npar,(char **)0), npar),
	mS(_mx2), mNPoint(22), mHasResolution(true), mMassDependent(false)
{
	mShiftPropagator = new std::vector<std::vector<GPUDataStream<float2> *> >[mNPoint];
	for (int i = 0; i != mNPoint; ++i) {
		mShiftPropagator[i] = std::vector<std::vector<GPUDataStream<float2> *> >(mS.GetList()->GetNSets(), std::vector<GPUDataStream<float2> *>());
		for	(int j = 0; j != mS.GetList()->GetNSets(); ++j){
			for (int k = 0; k != mS.GetList()->GetNBlocks(j); ++k){
				mShiftPropagator[i][j].push_back(0);
			}
		}
	}

	mResolution.push_back(std::make_pair(-_mean1-5.188000e+00f*_sigma1, _frac*8.121848e-07f));
	mResolution.push_back(std::make_pair(-_mean1-3.936166e+00f*_sigma1, _frac*1.956717e-04f));
	mResolution.push_back(std::make_pair(-_mean1-2.865126e+00f*_sigma1, _frac*6.720288e-03f));
	mResolution.push_back(std::make_pair(-_mean1-1.876039e+00f*_sigma1, _frac*6.613882e-02f));
	mResolution.push_back(std::make_pair(-_mean1-9.288696e-01f*_sigma1, _frac*2.422404e-01f));
	mResolution.push_back(std::make_pair(-_mean1+0.000000e+00f*_sigma1, _frac*3.694082e-01f));
	mResolution.push_back(std::make_pair(-_mean1+9.288696e-01f*_sigma1, _frac*2.422404e-01f));
	mResolution.push_back(std::make_pair(-_mean1+1.876039e+00f*_sigma1, _frac*6.613882e-02f));
	mResolution.push_back(std::make_pair(-_mean1+2.865126e+00f*_sigma1, _frac*6.720288e-03f));
	mResolution.push_back(std::make_pair(-_mean1+3.936166e+00f*_sigma1, _frac*1.956717e-04f));
	mResolution.push_back(std::make_pair(-_mean1+5.188000e+00f*_sigma1, _frac*8.121848e-07f));
	
	mResolution.push_back(std::make_pair(-_mean2-5.188000e+00f*_sigma2, (1-_frac)*8.121848e-07f));
	mResolution.push_back(std::make_pair(-_mean2-3.936166e+00f*_sigma2, (1-_frac)*1.956717e-04f));
	mResolution.push_back(std::make_pair(-_mean2-2.865126e+00f*_sigma2, (1-_frac)*6.720288e-03f));
	mResolution.push_back(std::make_pair(-_mean2-1.876039e+00f*_sigma2, (1-_frac)*6.613882e-02f));
	mResolution.push_back(std::make_pair(-_mean2-9.288696e-01f*_sigma2, (1-_frac)*2.422404e-01f));
	mResolution.push_back(std::make_pair(-_mean2+0.000000e+00f*_sigma2, (1-_frac)*3.694082e-01f));
	mResolution.push_back(std::make_pair(-_mean2+9.288696e-01f*_sigma2, (1-_frac)*2.422404e-01f));
	mResolution.push_back(std::make_pair(-_mean2+1.876039e+00f*_sigma2, (1-_frac)*6.613882e-02f));
	mResolution.push_back(std::make_pair(-_mean2+2.865126e+00f*_sigma2, (1-_frac)*6.720288e-03f));
	mResolution.push_back(std::make_pair(-_mean2+3.936166e+00f*_sigma2, (1-_frac)*1.956717e-04f));
	mResolution.push_back(std::make_pair(-_mean2+5.188000e+00f*_sigma2, (1-_frac)*8.121848e-07f));
}

// user-defined mass-dependent gaussian resolution
//GPUPropagatorResolution::GPUPropagatorResolution(GPUStreamScalar & _mx2, int _np, float * _resPars):
//	mS(_mx2), mNPoint(11), mHasResolution(true), mMassDependent(true)
//{
//	mResolution.push_back(std::make_pair(-5.188000e+00f, 8.121848e-07f));
//	mResolution.push_back(std::make_pair(-3.936166e+00f, 1.956717e-04f));
//	mResolution.push_back(std::make_pair(-2.865126e+00f, 6.720288e-03f));
//	mResolution.push_back(std::make_pair(-1.876039e+00f, 6.613882e-02f));
//	mResolution.push_back(std::make_pair(-9.288696e-01f, 2.422404e-01f));
//	mResolution.push_back(std::make_pair(0.000000e+00f, 3.694082e-01f));
//	mResolution.push_back(std::make_pair(9.288696e-01f, 2.422404e-01f));
//	mResolution.push_back(std::make_pair(1.876039e+00f, 6.613882e-02f));
//	mResolution.push_back(std::make_pair(2.865126e+00f, 6.720288e-03f));
//	mResolution.push_back(std::make_pair(3.936166e+00f, 1.956717e-04f));
//	mResolution.push_back(std::make_pair(5.188000e+00f, 8.121848e-07f));
//}

GPUPropagatorResolution::~GPUPropagatorResolution(void)
{
	//delete mShiftPropagator;
}
		
float GPUPropagatorResolution::GetShift(int num)
{
	assert(num >= 0 && num < mNPoint);
	if (mMassDependent) {
		return 0;
	}
	else
		return mResolution.at(num).first;
}

float GPUPropagatorResolution::GetWeight(int num)
{
	assert(num >= 0 && num < mNPoint);
	if (mMassDependent) {
		return 0;
	}
	else
		return mResolution.at(num).second;
}

void GPUPropagatorResolution::ResetResolution(unsigned int index)
{
	for (int r = 0; r != mNPoint; ++r) {
		if (index < (int) mShiftPropagator[r].size()) {
			for (int i = 0; i < (int) mShiftPropagator[r][index].size(); i++) {
				if (mShiftPropagator[r][index][i])
					delete mShiftPropagator[r][index][i];
				mShiftPropagator[r][index][i] = 0;
			}
		}
	}
}
	
void GPUPropagatorResolution::ResetResolution(unsigned int index, unsigned int block)
{
	for (int r = 0; r != mNPoint; ++r) {
		if (index < (int) mShiftPropagator[r].size()) {
			if (mShiftPropagator[r][index][block])
				delete mShiftPropagator[r][index][block];
			mShiftPropagator[r][index][block] = 0;
		}
	}
}

GPUDataStream<float2> * GPUPropagatorResolution::SumResolution(unsigned int index, unsigned int block)
{
	GPUDataStream<float2> * sum = new GPUDataStream<float2>(mS.GetList()->GetDeviceInterface(), mS.GetLength(index, block));
//        cout<<"GetNPoint()  "<<GetNPoint()<<endl;
	for (int i = 0; i != GetNPoint(); ++i) {

  //      cout<<"GetWeight(i)  "<<GetWeight(i)<<endl;
		float tag = (i==0 ? -1 : 1);
		kernelsumresolution(mS.GetList()->GetDeviceInterface(), tag, 
				GetShiftPropagator(i, index, block), 
				GetWeight(i),
				sum);
	}
	return sum;
}

GPUDataStream<float2> * GPUPropagatorResolution::ContractResolution(GPUPropagatorResolution & prop2, unsigned int index, unsigned int block)
{
	assert(IdenticalMx2(prop2));
	GPUDataStream<float2> * contract = new GPUDataStream<float2>(mS.GetList()->GetDeviceInterface(), mS.GetLength(index, block));
	for (int i = 0; i != GetNPoint(); ++i) {
//        cout<<"GetWeight(i)  "<<GetWeight(i)<<endl;

		float tag = (i==0 ? -1 : 1);
		kernelcontractresolution(mS.GetList()->GetDeviceInterface(), tag, 
				GetShiftPropagator(i, index, block), 
				prop2.GetShiftPropagator(i, index, block), 
				GetWeight(i),
				contract);
	}
	return contract;
}
