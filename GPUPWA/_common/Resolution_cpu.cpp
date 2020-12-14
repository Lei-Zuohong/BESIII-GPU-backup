/* Resolution_cpu.h  -- source for interface functions to the kernels
 defined in /besfs/groups/gpugroup/mintj/tmp/gpupwa/GPUPWA/Resolution.cl, the corresponding binaries are in binfiles/Resolution_cpu.bin
 THIS IS AN AUTO-GENERATED FILE, DO NOT EDIT */

#include"Resolution_cpu.h"
#include"KernelFile.h"
#include<vector>
#include<iostream>
#include<cassert>


namespace Resolution_CPU{

	 cl::Kernel * k_kernelsumresolution = 0;
	 cl::Kernel * k_kernelcontractresolution = 0;
	 cl::Kernel * k_kernelcontract_reswave_10 = 0;
	 cl::Kernel * k_kernelcontract_reswave_11 = 0;
	 cl::Kernel * k_kernelcontract_reswave_12 = 0;
	 cl::Kernel * k_kernelcontract_reswave_20 = 0;
	 cl::Kernel * k_kernelcontract_reswave_21 = 0;
	 cl::Kernel * k_kernelcontract_reswave_22 = 0;
	 cl::Kernel * k_kernelradcontract_reswave_10 = 0;
	 cl::Kernel * k_kernelradcontract_reswave_11 = 0;
	 cl::Kernel * k_kernelradcontract_reswave_12 = 0;
	 cl::Kernel * k_kernelradcontract_reswave_20 = 0;
	 cl::Kernel * k_kernelradcontract_reswave_21 = 0;
	 cl::Kernel * k_kernelradcontract_reswave_22 = 0;
	 int prepare_kernels(const DeviceInterface * dev){
		 FILE * input = fopen("binfiles/Resolution_cpu.bin","rb");
		 if(!input){
			 	std::cerr << "Loading binary file binfiles/Resolution_cpu.bin failed \n";
			 	return FAILURE;
		 }
		 fseek(input, 0L, SEEK_END);
		 size_t size = ftell(input);
		 rewind(input);
		 cl_int err;
		 char * binary = new char[size];
		 fread(binary, sizeof(char), size, input);
		 fclose(input);
		 cl::Program::Binaries binaries;
		 std::vector<int> binstatus;
		 for(unsigned int d=0; d < (*(dev->GetDevices())).size(); d++){
			 binstatus.push_back(0);
			 binaries.push_back(std::make_pair(binary, size));
		 }
		 cl::Program * bProgram = new cl::Program(*(dev->GetContext()),
		 						   *(dev->GetDevices()),
		 						   binaries,
		 						   &binstatus,
		 						   &err);
		 if (err != CL_SUCCESS) {
			 	std::cerr << "Loading binaries failed (" << err << ")\n";
			 	return FAILURE;
		 }
		 err = bProgram->build(*(dev->GetDevices()));
		 if (err != CL_SUCCESS) {
		 	std::cerr << "Program::build() failed (" << err << ")\n";
		 	return FAILURE;
		 }

		 k_kernelsumresolution = new cl::Kernel(*bProgram, "kernelsumresolution", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelsumresolution\n";
			 return FAILURE;
		 }

		 k_kernelcontractresolution = new cl::Kernel(*bProgram, "kernelcontractresolution", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelcontractresolution\n";
			 return FAILURE;
		 }

		 k_kernelcontract_reswave_10 = new cl::Kernel(*bProgram, "kernelcontract_reswave_10", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelcontract_reswave_10\n";
			 return FAILURE;
		 }

		 k_kernelcontract_reswave_11 = new cl::Kernel(*bProgram, "kernelcontract_reswave_11", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelcontract_reswave_11\n";
			 return FAILURE;
		 }

		 k_kernelcontract_reswave_12 = new cl::Kernel(*bProgram, "kernelcontract_reswave_12", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelcontract_reswave_12\n";
			 return FAILURE;
		 }

		 k_kernelcontract_reswave_20 = new cl::Kernel(*bProgram, "kernelcontract_reswave_20", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelcontract_reswave_20\n";
			 return FAILURE;
		 }

		 k_kernelcontract_reswave_21 = new cl::Kernel(*bProgram, "kernelcontract_reswave_21", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelcontract_reswave_21\n";
			 return FAILURE;
		 }

		 k_kernelcontract_reswave_22 = new cl::Kernel(*bProgram, "kernelcontract_reswave_22", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelcontract_reswave_22\n";
			 return FAILURE;
		 }

		 k_kernelradcontract_reswave_10 = new cl::Kernel(*bProgram, "kernelradcontract_reswave_10", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelradcontract_reswave_10\n";
			 return FAILURE;
		 }

		 k_kernelradcontract_reswave_11 = new cl::Kernel(*bProgram, "kernelradcontract_reswave_11", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelradcontract_reswave_11\n";
			 return FAILURE;
		 }

		 k_kernelradcontract_reswave_12 = new cl::Kernel(*bProgram, "kernelradcontract_reswave_12", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelradcontract_reswave_12\n";
			 return FAILURE;
		 }

		 k_kernelradcontract_reswave_20 = new cl::Kernel(*bProgram, "kernelradcontract_reswave_20", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelradcontract_reswave_20\n";
			 return FAILURE;
		 }

		 k_kernelradcontract_reswave_21 = new cl::Kernel(*bProgram, "kernelradcontract_reswave_21", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelradcontract_reswave_21\n";
			 return FAILURE;
		 }

		 k_kernelradcontract_reswave_22 = new cl::Kernel(*bProgram, "kernelradcontract_reswave_22", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelradcontract_reswave_22\n";
			 return FAILURE;
		 }

	 return 0;
	 }



	 int kernelsumresolution(const DeviceInterface * dev, float tag, Stream<float2> * prop, float weight, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 err = k_kernelsumresolution->setArg(0, tag);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((prop->GetEvent()));
		 err = k_kernelsumresolution->setArg(1, *(*prop)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelsumresolution->setArg(2, weight);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelsumresolution->setArg(3, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelsumresolution, cl::NullRange, cl::NDRange(prop->GetD1(),prop->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelcontractresolution(const DeviceInterface * dev, float tag, Stream<float2> * prop1, Stream<float2> * prop2, float weight, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 err = k_kernelcontractresolution->setArg(0, tag);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((prop1->GetEvent()));
		 err = k_kernelcontractresolution->setArg(1, *(*prop1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(prop2 != prop1)
			 evvec->push_back((prop2->GetEvent()));
		 err = k_kernelcontractresolution->setArg(2, *(*prop2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelcontractresolution->setArg(3, weight);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelcontractresolution->setArg(4, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelcontractresolution, cl::NullRange, cl::NDRange(prop1->GetD1(),prop1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelcontract_reswave_10(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float4> * tb, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelcontract_reswave_10->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelcontract_reswave_10->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelcontract_reswave_10->setArg(2, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelcontract_reswave_10->setArg(3, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelcontract_reswave_10, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelcontract_reswave_11(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float4> * tb, Stream<float2> * q1, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelcontract_reswave_11->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelcontract_reswave_11->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelcontract_reswave_11->setArg(2, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q1 != p1)
			 evvec->push_back((q1->GetEvent()));
		 err = k_kernelcontract_reswave_11->setArg(3, *(*q1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelcontract_reswave_11->setArg(4, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelcontract_reswave_11, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelcontract_reswave_12(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float4> * tb, Stream<float2> * q1, Stream<float2> * q2, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelcontract_reswave_12->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelcontract_reswave_12->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelcontract_reswave_12->setArg(2, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q1 != p1)
			 evvec->push_back((q1->GetEvent()));
		 err = k_kernelcontract_reswave_12->setArg(3, *(*q1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q2 != p1)
			if(q2 != q1)
				 evvec->push_back((q2->GetEvent()));
		 err = k_kernelcontract_reswave_12->setArg(4, *(*q2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelcontract_reswave_12->setArg(5, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelcontract_reswave_12, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelcontract_reswave_20(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float4> * tb, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelcontract_reswave_20->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelcontract_reswave_20->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(p2 != p1)
			 evvec->push_back((p2->GetEvent()));
		 err = k_kernelcontract_reswave_20->setArg(2, *(*p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelcontract_reswave_20->setArg(3, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelcontract_reswave_20->setArg(4, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelcontract_reswave_20, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelcontract_reswave_21(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float4> * tb, Stream<float2> * q1, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelcontract_reswave_21->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelcontract_reswave_21->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(p2 != p1)
			 evvec->push_back((p2->GetEvent()));
		 err = k_kernelcontract_reswave_21->setArg(2, *(*p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelcontract_reswave_21->setArg(3, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q1 != p1)
			if(q1 != p2)
				 evvec->push_back((q1->GetEvent()));
		 err = k_kernelcontract_reswave_21->setArg(4, *(*q1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelcontract_reswave_21->setArg(5, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelcontract_reswave_21, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelcontract_reswave_22(const DeviceInterface * dev, Stream<float4> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float4> * tb, Stream<float2> * q1, Stream<float2> * q2, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelcontract_reswave_22->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelcontract_reswave_22->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(p2 != p1)
			 evvec->push_back((p2->GetEvent()));
		 err = k_kernelcontract_reswave_22->setArg(2, *(*p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelcontract_reswave_22->setArg(3, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q1 != p1)
			if(q1 != p2)
				 evvec->push_back((q1->GetEvent()));
		 err = k_kernelcontract_reswave_22->setArg(4, *(*q1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q2 != p1)
			if(q2 != p2)
				if(q2 != q1)
					 evvec->push_back((q2->GetEvent()));
		 err = k_kernelcontract_reswave_22->setArg(5, *(*q2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelcontract_reswave_22->setArg(6, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelcontract_reswave_22, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelradcontract_reswave_10(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float44> * tb, Stream<float44> * gpp, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelradcontract_reswave_10->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelradcontract_reswave_10->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelradcontract_reswave_10->setArg(2, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(gpp != ta)
			if(gpp != tb)
				 evvec->push_back((gpp->GetEvent()));
		 err = k_kernelradcontract_reswave_10->setArg(3, *(*gpp)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelradcontract_reswave_10->setArg(4, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelradcontract_reswave_10, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelradcontract_reswave_11(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float44> * tb, Stream<float2> * q1, Stream<float44> * gpp, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelradcontract_reswave_11->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelradcontract_reswave_11->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelradcontract_reswave_11->setArg(2, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q1 != p1)
			 evvec->push_back((q1->GetEvent()));
		 err = k_kernelradcontract_reswave_11->setArg(3, *(*q1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(gpp != ta)
			if(gpp != tb)
				 evvec->push_back((gpp->GetEvent()));
		 err = k_kernelradcontract_reswave_11->setArg(4, *(*gpp)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelradcontract_reswave_11->setArg(5, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelradcontract_reswave_11, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelradcontract_reswave_12(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float44> * tb, Stream<float2> * q1, Stream<float2> * q2, Stream<float44> * gpp, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelradcontract_reswave_12->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelradcontract_reswave_12->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelradcontract_reswave_12->setArg(2, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q1 != p1)
			 evvec->push_back((q1->GetEvent()));
		 err = k_kernelradcontract_reswave_12->setArg(3, *(*q1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q2 != p1)
			if(q2 != q1)
				 evvec->push_back((q2->GetEvent()));
		 err = k_kernelradcontract_reswave_12->setArg(4, *(*q2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(gpp != ta)
			if(gpp != tb)
				 evvec->push_back((gpp->GetEvent()));
		 err = k_kernelradcontract_reswave_12->setArg(5, *(*gpp)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelradcontract_reswave_12->setArg(6, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelradcontract_reswave_12, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelradcontract_reswave_20(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float44> * tb, Stream<float44> * gpp, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelradcontract_reswave_20->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelradcontract_reswave_20->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(p2 != p1)
			 evvec->push_back((p2->GetEvent()));
		 err = k_kernelradcontract_reswave_20->setArg(2, *(*p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelradcontract_reswave_20->setArg(3, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(gpp != ta)
			if(gpp != tb)
				 evvec->push_back((gpp->GetEvent()));
		 err = k_kernelradcontract_reswave_20->setArg(4, *(*gpp)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelradcontract_reswave_20->setArg(5, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelradcontract_reswave_20, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelradcontract_reswave_21(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float44> * tb, Stream<float2> * q1, Stream<float44> * gpp, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelradcontract_reswave_21->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelradcontract_reswave_21->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(p2 != p1)
			 evvec->push_back((p2->GetEvent()));
		 err = k_kernelradcontract_reswave_21->setArg(2, *(*p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelradcontract_reswave_21->setArg(3, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q1 != p1)
			if(q1 != p2)
				 evvec->push_back((q1->GetEvent()));
		 err = k_kernelradcontract_reswave_21->setArg(4, *(*q1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(gpp != ta)
			if(gpp != tb)
				 evvec->push_back((gpp->GetEvent()));
		 err = k_kernelradcontract_reswave_21->setArg(5, *(*gpp)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelradcontract_reswave_21->setArg(6, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelradcontract_reswave_21, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

	 int kernelradcontract_reswave_22(const DeviceInterface * dev, Stream<float44> * ta, Stream<float2> * p1, Stream<float2> * p2, Stream<float44> * tb, Stream<float2> * q1, Stream<float2> * q2, Stream<float44> * gpp, Stream<float2> * output){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((ta->GetEvent()));
		 err = k_kernelradcontract_reswave_22->setArg(0, *(*ta)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 evvec->push_back((p1->GetEvent()));
		 err = k_kernelradcontract_reswave_22->setArg(1, *(*p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(p2 != p1)
			 evvec->push_back((p2->GetEvent()));
		 err = k_kernelradcontract_reswave_22->setArg(2, *(*p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(tb != ta)
			 evvec->push_back((tb->GetEvent()));
		 err = k_kernelradcontract_reswave_22->setArg(3, *(*tb)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q1 != p1)
			if(q1 != p2)
				 evvec->push_back((q1->GetEvent()));
		 err = k_kernelradcontract_reswave_22->setArg(4, *(*q1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(q2 != p1)
			if(q2 != p2)
				if(q2 != q1)
					 evvec->push_back((q2->GetEvent()));
		 err = k_kernelradcontract_reswave_22->setArg(5, *(*q2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(gpp != ta)
			if(gpp != tb)
				 evvec->push_back((gpp->GetEvent()));
		 err = k_kernelradcontract_reswave_22->setArg(6, *(*gpp)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelradcontract_reswave_22->setArg(7, *(*output)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelradcontract_reswave_22, cl::NullRange, cl::NDRange(ta->GetD1(),ta->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 output->SetEvent(event);

	 return 0;
	 }

}
