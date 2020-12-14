/* UserAmplitude_cpu.h  -- header for interface functions to the kernels
 defined in /besfs/groups/gpugroup/mintj/tmp/gpupwa/GPUPWA/UserAmplitude.cl, the corresponding binaries are in binfiles/UserAmplitude_cpu.bin
 THIS IS AN AUTO-GENERATED FILE, DO NOT EDIT */

#ifndef HHUserAmplitude_CPUHH
#define HHUserAmplitude_CPUHH


#include"Stream.h"
#include"DeviceInterface.h"


namespace UserAmplitude_CPU{

	 int prepare_kernels(const DeviceInterface * dev);

	 extern cl::Kernel * k_kernelusercontract;
	 extern cl::Kernel * k_kerneluserradiativecontract;
	 extern cl::Kernel * k_kerneluserpartialwave_trivial;
	 extern cl::Kernel * k_kerneluserradiativepartialwave_trivial;
	 extern cl::Kernel * k_kernel_gkk_gf0;
	 extern cl::Kernel * k_kernel_gkk_gf20;
	 extern cl::Kernel * k_kernel_gkk_gf21;
	 extern cl::Kernel * k_kernel_gkk_gf22;
	 extern cl::Kernel * k_kernel_kkpi_phi1pi;
	 extern cl::Kernel * k_kernel_kkpi_phi3pi;
	 int kernelusercontract(const DeviceInterface * dev, Stream<float4> * amp1, Stream<float4> * amp2, Stream<float2> * result);

	 int kerneluserradiativecontract(const DeviceInterface * dev, Stream<float44> * amp1, Stream<float44> * amp2, Stream<float44> * gpp, Stream<float2> * result);

	 int kerneluserpartialwave_trivial(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float4> * result);

	 int kerneluserradiativepartialwave_trivial(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float44> * result);

	 int kernel_gkk_gf0(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float44> * result);

	 int kernel_gkk_gf20(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float44> * result);

	 int kernel_gkk_gf21(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float44> * result);

	 int kernel_gkk_gf22(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float44> * result);

	 int kernel_kkpi_phi1pi(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float4> * result);

	 int kernel_kkpi_phi3pi(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float4> * result);

}
#endif
