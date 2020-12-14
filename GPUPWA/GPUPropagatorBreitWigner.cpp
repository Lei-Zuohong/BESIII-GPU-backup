#include "GPUPropagatorBreitWigner.h"
#ifdef USECPU
#include "Propagators_cpu.h"
using namespace Propagators_CPU;
#else
#include "Propagators.h"
using namespace Propagators_GPU;
#endif

#include <iostream>


GPUPropagatorBreitWigner::GPUPropagatorBreitWigner(char * name, GPUStreamScalar & _mx2):
                GPUPropagatorResolution(name, _mx2, std::vector<char**>(2,(char **)0), 2)
{       
        unsigned int sl = strlen(name);
        char ** p1 = new char*;
        char ** p2 = new char*;
        
        *p1 = new char[sl+8];
        *p2 = new char[sl+8];
        sprintf(*p1,"%s_%s",name,"mass");
        sprintf(*p2,"%s_%s",name,"width");
        mparnames[0] = p1;
        mparnames[1] = p2;
}

GPUPropagatorBreitWigner::GPUPropagatorBreitWigner(char * name, GPUStreamScalar & _mx2, float _mean, float _sigma):
                GPUPropagatorResolution(name, _mx2, std::vector<char**>(2,(char **)0), 2, _mean, _sigma)
{
        unsigned int sl = strlen(name);
        char ** p1 = new char*;
        char ** p2 = new char*;

        *p1 = new char[sl+8];
        *p2 = new char[sl+8];
        sprintf(*p1,"%s_%s",name,"mass");
        sprintf(*p2,"%s_%s",name,"width");
        mparnames[0] = p1;
        mparnames[1] = p2;
}




GPUPropagatorBreitWigner::GPUPropagatorBreitWigner(char * name, GPUStreamScalar & _mx2, float _mean1, float _sigma1, float _mean2, float _sigma2, float _frac):
                GPUPropagatorResolution(name, _mx2, std::vector<char**>(2,(char **)0), 2, _mean1, _sigma1, _mean2, _sigma2, _frac)
{
        unsigned int sl = strlen(name);
        char ** p1 = new char*;
        char ** p2 = new char*;

        *p1 = new char[sl+8];
        *p2 = new char[sl+8];
        sprintf(*p1,"%s_%s",name,"mass");
        sprintf(*p2,"%s_%s",name,"width");
        mparnames[0] = p1;
        mparnames[1] = p2;
}




/*
GPUPropagatorBreitWigner::GPUPropagatorBreitWigner(char * name,GPUStreamScalar & _mx2):
		GPUBasicPropagator(name, _mx2, std::vector<char**>(2,(char **)0), 2)
{
	unsigned int sl = strlen(name);
	char ** p1 = new char*;
	char ** p2 = new char*;

	*p1 = new char[sl+8];
	*p2 = new char[sl+8];
	sprintf(*p1,"%s_%s",name,"mass");
	sprintf(*p2,"%s_%s",name,"width");
	mparnames[0] = p1;
	mparnames[1] = p2;
}
*/
GPUPropagatorBreitWigner::~GPUPropagatorBreitWigner(void)
{
}


GPUDataStream<float2> * GPUPropagatorBreitWigner::operator()(int index, int block){
	assert(index < mx2.GetNStreams());
	assert((int)mstream.size() > index);
	if(CacheValid(index)){
		//cout << "Cache valid!" << endl;
		if(mstream[index][block]){
			return mstream[index][block];
		}
	} else {
		InvalidateCache(index);
	}
	Stream<float> * mx2stream = mx2(index, block);

	mstream[index][block] = new GPUDataStream<float2>(mList->GetDeviceInterface(), mx2.GetLength(index, block));
	double mass = GetParameterValue(0);
	double width = GetParameterValue(1);

	//cout << index << " : " << block << " Calling with mass: " << mass <<" and width " << width << endl;

	//kernelbreitwignerRes(mList->GetDeviceInterface(), mx2stream,0,mass, mass*mass, width, mstream[index][block]);
	kernelmassdependentbreitwigner1Reskk(mList->GetDeviceInterface(), mx2stream,0,mass, mass*mass, width, mstream[index][block]);
	//kernelbreitwigner(mList->GetDeviceInterface(), mx2stream,mass, mass*mass, width, mstream[index][block]);
	SetLastValue(index, 0, mass);
	SetLastValue(index, 1, width);
	return mstream[index][block];

};


/// Get a stream of the derivative of the propagator with regards to the mass
GPUDataStream<float2> * GPUPropagatorBreitWigner::GetMassDerivative(int index, int block){
	assert(index < mx2.GetNStreams());
	assert((int)mDmass.size() > index);
	if(CacheValid(index)){
		if(mDmass[index][block])
			return mDmass[index][block];
	} else {
		InvalidateCache(index);
	}
	Stream<float> * mx2stream = mx2(index,block);

	mDmass[index][block] = new GPUDataStream<float2>(mList->GetDeviceInterface(), mx2.GetLength(index, block));
	double mass = GetParameterValue(0);
	double width = GetParameterValue(1);

	kerneldbreitwignerdmass(mList->GetDeviceInterface(), mx2stream,mass, mass*mass, width, mDmass[index][block]);
	SetLastValue(index, 0, mass);
	SetLastValue(index, 1, width);
	return mDmass[index][block];

}
/// Get a stream of the derivative of the propagator with regards to the width
GPUDataStream<float2> * GPUPropagatorBreitWigner::GetWidthDerivative(int index, int block){
	assert(index < mx2.GetNStreams());
	assert((int)mDwidth.size() > index);
	if(CacheValid(index)){
		if(mDwidth[index][block])
			return mDwidth[index][block];
	} else {
		InvalidateCache(index);
	}
	Stream<float> * mx2stream = mx2(index, block);

	mDwidth[index][block] = new GPUDataStream<float2>(mList->GetDeviceInterface(), mx2.GetLength(index, block));
	double mass = GetParameterValue(0);
	double width = GetParameterValue(1);
	kerneldbreitwignerdwidth(mList->GetDeviceInterface(), mx2stream,mass, mass*mass, width, mDwidth[index][block]);
	SetLastValue(index, 0, mass);
	SetLastValue(index, 1, width);
	return mDwidth[index][block];
}

// Get a stream of the derivative of the |propagator|^2 with regards to the mass
GPUDataStream<float> * GPUPropagatorBreitWigner::GetAbsMassDerivative(int index, int block){
	assert(index < mx2.GetNStreams());
	assert((int)mAbsDmass.size() > index);
	if(CacheValid(index)){
		if(mAbsDmass[index][block])
			return mAbsDmass[index][block];
	} else {
		InvalidateCache(index);
	}
	Stream<float> * mx2stream = mx2(index, block);

	mAbsDmass[index][block] = new GPUDataStream<float>(mList->GetDeviceInterface(), mx2.GetLength(index, block));
	double mass = GetParameterValue(0);
	double width = GetParameterValue(1);
	kerneldbw2dmass(mList->GetDeviceInterface(), mx2stream,mass, mass*mass, width, mAbsDmass[index][block]);
	SetLastValue(index, 0, mass);
	SetLastValue(index, 1, width);
	return mAbsDmass[index][block];

}
/// Get a stream of the derivative of the |propagator|^2 with regards to the width
GPUDataStream<float> * GPUPropagatorBreitWigner::GetAbsWidthDerivative(int index, int block){
	assert(index < mx2.GetNStreams());
	assert((int)mAbsDwidth.size() > index);
	if(CacheValid(index)){
		if(mAbsDwidth[index][block])
			return mAbsDwidth[index][block];
	} else {
		InvalidateCache(index);
	}
	Stream<float> * mx2stream = mx2(index, block);

	mAbsDwidth[index][block] = new GPUDataStream<float>(mList->GetDeviceInterface(), mx2.GetLength(index, block));
	double mass = GetParameterValue(0);
	double width = GetParameterValue(1);
	kerneldbw2dwidth(mList->GetDeviceInterface(), mx2stream,mass, mass*mass, width, mAbsDwidth[index][block]);
	SetLastValue(index, 0, mass);
	SetLastValue(index, 1, width);
	return mAbsDwidth[index][block];

}

GPUDataStream<float2> * GPUPropagatorBreitWigner::GetShiftPropagator(unsigned int num, unsigned int index, unsigned int block)
{
//        cout<<"num  "<<num<<"   mNPoint   "<<mNPoint<<endl;
        assert(num >= 0 && num < mNPoint);
        assert(index < mx2.GetNStreams());
        assert((int)mShiftPropagator[num].size() > index);

        if (CacheValid(index)) {
                if (mShiftPropagator[num][index][block]) {
                        return mShiftPropagator[num][index][block];
                }
        }
        else {
                InvalidateCache(index);
        }


        float shift = GetShift(num);
  //      cout<<"shift  "<<shift<<endl;
        Stream <float> * mx2stream = mx2(index, block);
        mShiftPropagator[num][index][block] = new GPUDataStream<float2> (mList->GetDeviceInterface(), mx2.GetLength(index, block));
        double mass = GetParameterValue(0);
        double width = GetParameterValue(1);
        //kernelbreitwignerRes(mList->GetDeviceInterface(), mx2stream, shift, mass, mass*mass, width, mShiftPropagator[num][index][block]);
        kernelmassdependentbreitwigner1Reskk(mList->GetDeviceInterface(), mx2stream, shift, mass, mass*mass, width, mShiftPropagator[num][index][block]);
        SetLastValue(index, 0, mass);
        SetLastValue(index, 1, width);
        return mShiftPropagator[num][index][block];
}
