/* Resolution_cpu.h  -- header for interface functions to the kernels
 defined in /besfs/groups/gpugroup/mintj/tmp/gpupwa/GPUPWA/Resolution.cl, the corresponding binaries are in binfiles/Resolution_cpu.bin
 THIS IS AN AUTO-GENERATED FILE, DO NOT EDIT */

#ifndef HHResolution_CPUHH
#define HHResolution_CPUHH


#include"Stream.h"
#include"DeviceInterface.h"


namespace Resolution_CPU{

	 int prepare_kernels(const DeviceInterface * dev);

	 extern cl::Kernel * k_kernelsumresolution;
	 extern cl::Kernel * k_kernelcontractresolution;
	 extern cl::Kernel * k_kernelcontract_reswave_10;
	 extern cl::Kernel * k_kernelcontract_reswave_11;
	 extern cl::Kernel * k_kernelcontract_reswave_12;
	 extern cl::Kernel * k_kernelcontract_reswave_20;
	 extern cl::Kernel * k_kernelcontract_reswave_21;
	 extern cl::Kernel * k_kernelcontract_reswave_22;
	 extern cl::Kernel * k_kernelradcontract_reswave_10;
	 extern cl::Kernel * k_kernelradcontract_reswave_11;
	 extern cl::Kernel * k_kernelradcontract_reswave_12;
	 extern cl::Kernel * k_kernelradcontract_reswave_20;
	 extern cl::Kernel * k_kernelradcontract_reswave_21;
	 extern cl::Kernel * k_kernelradcontract_reswave_22;
	 int kernelsumresolution(const DeviceInterface * dev, float tag, Stream<float2> * prop, float weight, Stream<float2> * output);

	 int kernelcontractresolution(const DeviceInterface * dev, float tag, Stream<float2> * prop1, Stream<float2> * prop2, float weight, Stream<float2> * output);

	 int kernelcontract_reswave_10(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float4> * tb, Stream<float2> * output);

	 int kernelcontract_reswave_11(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float4> * tb, Stream<float2> * q1, Stream<float2> * output);

	 int kernelcontract_reswave_12(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float4> * tb, Stream<float2> * q1, Stream<float2> * q2, Stream<float2> * output);

	 int kernelcontract_reswave_20(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float4> * tb, Stream<float2> * output);

	 int kernelcontract_reswave_21(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float4> * tb, Stream<float2> * q1, Stream<float2> * output);

	 int kernelcontract_reswave_22(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float4> * tb, Stream<float2> * q1, Stream<float2> * q2, Stream<float2> * output);

	 int kernelradcontract_reswave_10(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float44> * tb, Stream<float44> * gpp, Stream<float2> * output);

	 int kernelradcontract_reswave_11(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float44> * tb, Stream<float2> * q1, Stream<float44> * gpp, Stream<float2> * output);

	 int kernelradcontract_reswave_12(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float44> * tb, Stream<float2> * q1, Stream<float2> * q2, Stream<float44> * gpp, Stream<float2> * output);

	 int kernelradcontract_reswave_20(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float44> * tb, Stream<float44> * gpp, Stream<float2> * output);

	 int kernelradcontract_reswave_21(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float44> * tb, Stream<float2> * q1, Stream<float44> * gpp, Stream<float2> * output);

	 int kernelradcontract_reswave_22(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float44> * tb, Stream<float2> * q1, Stream<float2> * q2, Stream<float44> * gpp, Stream<float2> * output);

}
#endif
