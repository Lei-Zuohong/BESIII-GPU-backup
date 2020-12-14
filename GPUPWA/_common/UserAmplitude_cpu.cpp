/* UserAmplitude_cpu.h  -- source for interface functions to the kernels
 defined in /besfs/groups/gpugroup/mintj/tmp/gpupwa/GPUPWA/UserAmplitude.cl, the corresponding binaries are in binfiles/UserAmplitude_cpu.bin
 THIS IS AN AUTO-GENERATED FILE, DO NOT EDIT */

#include"UserAmplitude_cpu.h"
#include"KernelFile.h"
#include<vector>
#include<iostream>
#include<cassert>


namespace UserAmplitude_CPU{

	 cl::Kernel * k_kernelusercontract = 0;
	 cl::Kernel * k_kerneluserradiativecontract = 0;
	 cl::Kernel * k_kerneluserpartialwave_trivial = 0;
	 cl::Kernel * k_kerneluserradiativepartialwave_trivial = 0;
	 cl::Kernel * k_kernel_gkk_gf0 = 0;
	 cl::Kernel * k_kernel_gkk_gf20 = 0;
	 cl::Kernel * k_kernel_gkk_gf21 = 0;
	 cl::Kernel * k_kernel_gkk_gf22 = 0;
	 cl::Kernel * k_kernel_kkpi_phi1pi = 0;
	 cl::Kernel * k_kernel_kkpi_phi3pi = 0;
	 int prepare_kernels(const DeviceInterface * dev){
		 FILE * input = fopen("binfiles/UserAmplitude_cpu.bin","rb");
		 if(!input){
			 	std::cerr << "Loading binary file binfiles/UserAmplitude_cpu.bin failed \n";
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

		 k_kernelusercontract = new cl::Kernel(*bProgram, "kernelusercontract", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernelusercontract\n";
			 return FAILURE;
		 }

		 k_kerneluserradiativecontract = new cl::Kernel(*bProgram, "kerneluserradiativecontract", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kerneluserradiativecontract\n";
			 return FAILURE;
		 }

		 k_kerneluserpartialwave_trivial = new cl::Kernel(*bProgram, "kerneluserpartialwave_trivial", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kerneluserpartialwave_trivial\n";
			 return FAILURE;
		 }

		 k_kerneluserradiativepartialwave_trivial = new cl::Kernel(*bProgram, "kerneluserradiativepartialwave_trivial", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kerneluserradiativepartialwave_trivial\n";
			 return FAILURE;
		 }

		 k_kernel_gkk_gf0 = new cl::Kernel(*bProgram, "kernel_gkk_gf0", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernel_gkk_gf0\n";
			 return FAILURE;
		 }

		 k_kernel_gkk_gf20 = new cl::Kernel(*bProgram, "kernel_gkk_gf20", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernel_gkk_gf20\n";
			 return FAILURE;
		 }

		 k_kernel_gkk_gf21 = new cl::Kernel(*bProgram, "kernel_gkk_gf21", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernel_gkk_gf21\n";
			 return FAILURE;
		 }

		 k_kernel_gkk_gf22 = new cl::Kernel(*bProgram, "kernel_gkk_gf22", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernel_gkk_gf22\n";
			 return FAILURE;
		 }

		 k_kernel_kkpi_phi1pi = new cl::Kernel(*bProgram, "kernel_kkpi_phi1pi", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernel_kkpi_phi1pi\n";
			 return FAILURE;
		 }

		 k_kernel_kkpi_phi3pi = new cl::Kernel(*bProgram, "kernel_kkpi_phi3pi", &err);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel::Kernel() failed (" << err << ") for kernel kernel_kkpi_phi3pi\n";
			 return FAILURE;
		 }

	 return 0;
	 }



	 int kernelusercontract(const DeviceInterface * dev, Stream<float4> * amp1, Stream<float4> * amp2, Stream<float2> * result){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((amp1->GetEvent()));
		 err = k_kernelusercontract->setArg(0, *(*amp1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(amp2 != amp1)
			 evvec->push_back((amp2->GetEvent()));
		 err = k_kernelusercontract->setArg(1, *(*amp2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernelusercontract->setArg(2, *(*result)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernelusercontract, cl::NullRange, cl::NDRange(amp1->GetD1(),amp1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 result->SetEvent(event);

	 return 0;
	 }

	 int kerneluserradiativecontract(const DeviceInterface * dev, Stream<float44> * amp1, Stream<float44> * amp2, Stream<float44> * gpp, Stream<float2> * result){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((amp1->GetEvent()));
		 err = k_kerneluserradiativecontract->setArg(0, *(*amp1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(amp2 != amp1)
			 evvec->push_back((amp2->GetEvent()));
		 err = k_kerneluserradiativecontract->setArg(1, *(*amp2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(gpp != amp1)
			if(gpp != amp2)
				 evvec->push_back((gpp->GetEvent()));
		 err = k_kerneluserradiativecontract->setArg(2, *(*gpp)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kerneluserradiativecontract->setArg(3, *(*result)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kerneluserradiativecontract, cl::NullRange, cl::NDRange(amp1->GetD1(),amp1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 result->SetEvent(event);

	 return 0;
	 }

	 int kerneluserpartialwave_trivial(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float4> * result){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((_p1->GetEvent()));
		 err = k_kerneluserpartialwave_trivial->setArg(0, *(*_p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p2 != _p1)
			 evvec->push_back((_p2->GetEvent()));
		 err = k_kerneluserpartialwave_trivial->setArg(1, *(*_p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p3 != _p1)
			if(_p3 != _p2)
				 evvec->push_back((_p3->GetEvent()));
		 err = k_kerneluserpartialwave_trivial->setArg(2, *(*_p3)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kerneluserpartialwave_trivial->setArg(3, mass);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kerneluserpartialwave_trivial->setArg(4, width);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kerneluserpartialwave_trivial->setArg(5, *(*result)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kerneluserpartialwave_trivial, cl::NullRange, cl::NDRange(_p1->GetD1(),_p1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 result->SetEvent(event);

	 return 0;
	 }

	 int kerneluserradiativepartialwave_trivial(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float44> * result){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((_p1->GetEvent()));
		 err = k_kerneluserradiativepartialwave_trivial->setArg(0, *(*_p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p2 != _p1)
			 evvec->push_back((_p2->GetEvent()));
		 err = k_kerneluserradiativepartialwave_trivial->setArg(1, *(*_p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p3 != _p1)
			if(_p3 != _p2)
				 evvec->push_back((_p3->GetEvent()));
		 err = k_kerneluserradiativepartialwave_trivial->setArg(2, *(*_p3)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kerneluserradiativepartialwave_trivial->setArg(3, mass);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kerneluserradiativepartialwave_trivial->setArg(4, width);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kerneluserradiativepartialwave_trivial->setArg(5, *(*result)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kerneluserradiativepartialwave_trivial, cl::NullRange, cl::NDRange(_p1->GetD1(),_p1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 result->SetEvent(event);

	 return 0;
	 }

	 int kernel_gkk_gf0(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float44> * result){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((_p1->GetEvent()));
		 err = k_kernel_gkk_gf0->setArg(0, *(*_p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p2 != _p1)
			 evvec->push_back((_p2->GetEvent()));
		 err = k_kernel_gkk_gf0->setArg(1, *(*_p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p3 != _p1)
			if(_p3 != _p2)
				 evvec->push_back((_p3->GetEvent()));
		 err = k_kernel_gkk_gf0->setArg(2, *(*_p3)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf0->setArg(3, mass);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf0->setArg(4, width);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf0->setArg(5, *(*result)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernel_gkk_gf0, cl::NullRange, cl::NDRange(_p1->GetD1(),_p1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 result->SetEvent(event);

	 return 0;
	 }

	 int kernel_gkk_gf20(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float44> * result){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((_p1->GetEvent()));
		 err = k_kernel_gkk_gf20->setArg(0, *(*_p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p2 != _p1)
			 evvec->push_back((_p2->GetEvent()));
		 err = k_kernel_gkk_gf20->setArg(1, *(*_p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p3 != _p1)
			if(_p3 != _p2)
				 evvec->push_back((_p3->GetEvent()));
		 err = k_kernel_gkk_gf20->setArg(2, *(*_p3)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf20->setArg(3, mass);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf20->setArg(4, width);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf20->setArg(5, *(*result)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernel_gkk_gf20, cl::NullRange, cl::NDRange(_p1->GetD1(),_p1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 result->SetEvent(event);

	 return 0;
	 }

	 int kernel_gkk_gf21(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float44> * result){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((_p1->GetEvent()));
		 err = k_kernel_gkk_gf21->setArg(0, *(*_p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p2 != _p1)
			 evvec->push_back((_p2->GetEvent()));
		 err = k_kernel_gkk_gf21->setArg(1, *(*_p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p3 != _p1)
			if(_p3 != _p2)
				 evvec->push_back((_p3->GetEvent()));
		 err = k_kernel_gkk_gf21->setArg(2, *(*_p3)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf21->setArg(3, mass);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf21->setArg(4, width);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf21->setArg(5, *(*result)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernel_gkk_gf21, cl::NullRange, cl::NDRange(_p1->GetD1(),_p1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 result->SetEvent(event);

	 return 0;
	 }

	 int kernel_gkk_gf22(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float44> * result){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((_p1->GetEvent()));
		 err = k_kernel_gkk_gf22->setArg(0, *(*_p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p2 != _p1)
			 evvec->push_back((_p2->GetEvent()));
		 err = k_kernel_gkk_gf22->setArg(1, *(*_p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p3 != _p1)
			if(_p3 != _p2)
				 evvec->push_back((_p3->GetEvent()));
		 err = k_kernel_gkk_gf22->setArg(2, *(*_p3)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf22->setArg(3, mass);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf22->setArg(4, width);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_gkk_gf22->setArg(5, *(*result)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernel_gkk_gf22, cl::NullRange, cl::NDRange(_p1->GetD1(),_p1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 result->SetEvent(event);

	 return 0;
	 }

	 int kernel_kkpi_phi1pi(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float4> * result){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((_p1->GetEvent()));
		 err = k_kernel_kkpi_phi1pi->setArg(0, *(*_p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p2 != _p1)
			 evvec->push_back((_p2->GetEvent()));
		 err = k_kernel_kkpi_phi1pi->setArg(1, *(*_p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p3 != _p1)
			if(_p3 != _p2)
				 evvec->push_back((_p3->GetEvent()));
		 err = k_kernel_kkpi_phi1pi->setArg(2, *(*_p3)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_kkpi_phi1pi->setArg(3, mass);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_kkpi_phi1pi->setArg(4, width);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_kkpi_phi1pi->setArg(5, *(*result)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernel_kkpi_phi1pi, cl::NullRange, cl::NDRange(_p1->GetD1(),_p1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 result->SetEvent(event);

	 return 0;
	 }

	 int kernel_kkpi_phi3pi(const DeviceInterface * dev, Stream<float4> * _p1, Stream<float4> * _p2, Stream<float4> * _p3, float mass, float width, Stream<float4> * result){

		 cl_int err;
		 std::vector<cl::Event> *evvec = new std::vector<cl::Event>();
		 cl::Event event = cl::Event();

		 evvec->push_back((_p1->GetEvent()));
		 err = k_kernel_kkpi_phi3pi->setArg(0, *(*_p1)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p2 != _p1)
			 evvec->push_back((_p2->GetEvent()));
		 err = k_kernel_kkpi_phi3pi->setArg(1, *(*_p2)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		if(_p3 != _p1)
			if(_p3 != _p2)
				 evvec->push_back((_p3->GetEvent()));
		 err = k_kernel_kkpi_phi3pi->setArg(2, *(*_p3)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_kkpi_phi3pi->setArg(3, mass);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_kkpi_phi3pi->setArg(4, width);
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = k_kernel_kkpi_phi3pi->setArg(5, *(*result)());
		 if (err != CL_SUCCESS) {
			 std::cerr << "Kernel.SetArg() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 err = dev->GetQueue()->enqueueNDRangeKernel(*k_kernel_kkpi_phi3pi, cl::NullRange, cl::NDRange(_p1->GetD1(),_p1->GetD2()), cl::NullRange, evvec, &event );
		 if (err != CL_SUCCESS) {
			 std::cerr << "CommandQueue::enqueueNDRangeKernel() failed (" << err << ")\n";
			 assert(0);
			 return FAILURE;
		 }

		 result->SetEvent(event);

	 return 0;
	 }

}
