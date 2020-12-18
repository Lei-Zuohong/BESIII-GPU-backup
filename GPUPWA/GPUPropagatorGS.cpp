#include "GPUPropagatorGS.h"
#ifdef USECPU
#include "Propagators_cpu.h"
using namespace Propagators_CPU;
#else
#include "Propagators.h"
using namespace Propagators_GPU;
#endif

GPUPropagatorGS::GPUPropagatorGS(char *name,
								 GPUStreamScalar &_mx2,
								 float _m_1,
								 float _m_2,
								 float _r,
								 int _l) : GPUBasicPropagator(name,
															  _mx2,
															  std::vector<char **>(2, (char **)0), 2),
										   m1_2(_m_1 * _m_1),
										   m2_2(_m_2 * _m_2),
										   m_r(_r),
										   m_l(_l)
{
	unsigned int sl = strlen(name);
	char **p1 = new char *;
	char **p2 = new char *;
	*p1 = new char[sl + 8];
	*p2 = new char[sl + 8];
	sprintf(*p1, "%s_%s", name, "mass");
	sprintf(*p2, "%s_%s", name, "width");
	mparnames[0] = p1;
	mparnames[1] = p2;
}
GPUPropagatorGS::~GPUPropagatorGS(void)
{
}
GPUDataStream<float2> *GPUPropagatorGS::operator()(int index, int block)
{
	assert(index < mx2.GetNStreams());
	assert((int)mstream.size() > index);
	if (CacheValid(index))
	{
		if (mstream[index][block])
		{
			return mstream[index][block];
		}
	}
	else
	{
		InvalidateCache(index);
	}
	Stream<float> *mx2stream = mx2(index, block);

	mstream[index][block] = new GPUDataStream<float2>(mList->GetDeviceInterface(), mx2.GetLength(index, block));
	double mass = GetParameterValue(0);
	double width = GetParameterValue(1);

	kernelGS(mList->GetDeviceInterface(), mx2stream, mass, width, m1_2, m2_2, m_r, m_l, mstream[index][block]);
	SetLastValue(index, 0, mass);
	SetLastValue(index, 1, width);
	return mstream[index][block];
};
