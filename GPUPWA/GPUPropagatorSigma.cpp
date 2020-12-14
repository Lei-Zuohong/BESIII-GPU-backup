#include "GPUPropagatorSigma.h"
#include <cassert>
#include <cmath>
#ifdef USECPU
#include "Propagators_cpu.h"
using namespace Propagators_CPU;
#else
#include "Propagators.h"
using namespace Propagators_GPU;
#endif



GPUPropagatorSigma::GPUPropagatorSigma(char * name, GPUStreamScalar & _mx2):
		GPUBasicPropagator(name, _mx2, std::vector<char**>(2,(char **)0), 2){

        unsigned int sl = strlen(name);
        char ** p1 = new char*;
        char ** p2 = new char*;

        *p1 = new char[sl+8];
        *p2 = new char[sl+8];
        sprintf(*p1,"%s_%s",name,"mass");
        sprintf(*p2,"%s_%s",name,"gf");
        mparnames[0] = p1;
        mparnames[1] = p2;


}

GPUPropagatorSigma::~GPUPropagatorSigma(void)
{
}


GPUDataStream<float2> * GPUPropagatorSigma::operator()(int index, int block){
	assert(index < mx2.GetNStreams());
	assert((int)mstream.size() > index);
	if(CacheValid(index)){
		//cout << "Cache valid!" << endl;
		if(mstream[index][block]){
			return mstream[index][block];
		}
	} else {
		//cout << "Cache invalid!" << endl;
		InvalidateCache(index);
	}
	Stream<float> * mx2stream = mx2(index, block);

	mstream[index][block] = new GPUDataStream<float2>(mList->GetDeviceInterface(), mx2.GetLength(index, block));
        double mass = GetParameterValue(0);
        double width = GetParameterValue(1);



	//cout << index << " : " << block << " Calling with mass: " << mass <<" , g1: " << g1 << " and g2: " << g2 << endl;

	kernelsigma(mList->GetDeviceInterface(), mx2stream,mass,width,mstream[index][block]);
        SetLastValue(index, 0, mass);
        SetLastValue(index, 1, width);


	return mstream[index][block];
};
