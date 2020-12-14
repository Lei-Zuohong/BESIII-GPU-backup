/*/// \file Resolution.br*/
#pragma OPENCL EXTENSION cl_amd_fp64 : enable

#ifndef F44STRUCTS
#define F44STRUCTS

typedef struct {
  float4  c;
  float4  d;
  float4  e;
  float4  f;
} float44;

#endif

/* resolution */
__kernel void kernelsumresolution(float tag, __global float2 * prop, float weight, __global out float2 * output)
{
	uint i = get_global_id(0);
	if (tag < 0) {
		output[i].x = 0;
		output[i].y = 0;
	}
	
	output[i].x += weight*prop[i].x;
	output[i].y += weight*prop[i].y;
}

__kernel void kernelcontractresolution(float tag, __global float2 * prop1, __global float2 * prop2, float weight, __global out float2 * output)
{
	uint i = get_global_id(0);
	if (tag < 0) {
		output[i].x = 0;
		output[i].y = 0;
	}
	
	output[i].x += weight*(prop1[i].x*prop2[i].x + prop1[i].y*prop2[i].y);
	output[i].y += weight*(-prop1[i].x*prop2[i].y + prop1[i].y*prop2[i].x);
}

/* non-radiative */
__kernel void kernelcontract_reswave_10(__global float4 * ta, __global float2 * p1, __global float4 * tb, __global out float2 * output)
{
	uint i = get_global_id(0);
	float ts = 0.5f*(ta[i].x*tb[i].x + ta[i].y*tb[i].y);
	output[i].x = ts * (p1[i].x);
	output[i].y = ts * (p1[i].y);
}

__kernel void kernelcontract_reswave_11(__global float4 * ta, __global float2 * p1, __global float4 * tb, __global float2 * q1, __global out float2 * output)
{
	uint i = get_global_id(0);
	float ts = 0.5f*(ta[i].x*tb[i].x + ta[i].y*tb[i].y);
	output[i].x = ts * (p1[i].x*q1[i].x - p1[i].y*-q1[i].y);
	output[i].y = ts * (p1[i].x*-q1[i].y + p1[i].y*q1[i].x);
}

__kernel void kernelcontract_reswave_12(__global float4 * ta, __global float2 * p1, __global float4 * tb, __global float2 * q1, __global float2 * q2, __global out float2 * output)
{
	uint i = get_global_id(0);
	float ts = 0.5f*(ta[i].x*tb[i].x + ta[i].y*tb[i].y);
	output[i].x = ts * (p1[i].x*(q1[i].x*q2[i].x - q1[i].y*q2[i].y) - p1[i].y*-(q1[i].x*q2[i].y + q1[i].y*q2[i].x));
	output[i].y = ts * (p1[i].x*-(q1[i].x*q2[i].y + q1[i].y*q2[i].x) + p1[i].y*(q1[i].x*q2[i].x - q1[i].y*q2[i].y));
}

__kernel void kernelcontract_reswave_20(__global float4 * ta, __global float2 * p1, __global float2 * p2, __global float4 * tb, __global out float2 * output)
{
	uint i = get_global_id(0);
	float ts = 0.5f*(ta[i].x*tb[i].x + ta[i].y*tb[i].y);
	output[i].x = ts * (p1[i].x*p2[i].x - p1[i].y*p2[i].y);
	output[i].y = ts * (p1[i].x*p2[i].y + p1[i].y*p2[i].x);
}

__kernel void kernelcontract_reswave_21(__global float4 * ta, __global float2 * p1, __global float2 * p2, __global float4 * tb, __global float2 * q1, __global out float2 * output)
{
	uint i = get_global_id(0);
	float ts = 0.5f*(ta[i].x*tb[i].x + ta[i].y*tb[i].y);
	output[i].x = ts * ((p1[i].x*p2[i].x - p1[i].y*p2[i].y)*q1[i].x - (p1[i].x*p2[i].y + p1[i].y*p2[i].x)*-q1[i].y);
	output[i].y = ts * ((p1[i].x*p2[i].x - p1[i].y*p2[i].y)*-q1[i].y + (p1[i].x*p2[i].y + p1[i].y*p2[i].x)*q1[i].x);
}

__kernel void kernelcontract_reswave_22(__global float4 * ta, __global float2 * p1, __global float2 * p2, __global float4 * tb, __global float2 * q1, __global float2 * q2, __global out float2 * output)
{
	uint i = get_global_id(0);
	float ts = 0.5f*(ta[i].x*tb[i].x + ta[i].y*tb[i].y);
	output[i].x = ts * ((p1[i].x*p2[i].x - p1[i].y*p2[i].y)*(q1[i].x*q2[i].x - q1[i].y*q2[i].y) - (p1[i].x*p2[i].y + p1[i].y*p2[i].x)*-(q1[i].x*q2[i].y + q1[i].y*q2[i].x));
	output[i].y = ts * ((p1[i].x*p2[i].x - p1[i].y*p2[i].y)*-(q1[i].x*q2[i].y + q1[i].y*q2[i].x) + (p1[i].x*p2[i].y + p1[i].y*p2[i].x)*(q1[i].x*q2[i].x - q1[i].y*q2[i].y));
}

/* radiative */
__kernel void kernelradcontract_reswave_10(__global float44 * ta, __global float2 * p1, __global float44 * tb, __global float44 * gpp, __global out float2 * output)
{
	uint i = get_global_id(0);
	float44 sumoneindex;
	float4 sumtwoindex;
	sumoneindex.c.x = ta[i].c.x*gpp[i].c.x + ta[i].c.y*gpp[i].d.x + ta[i].c.z*gpp[i].e.x;
	sumoneindex.c.y = ta[i].c.x*gpp[i].c.y + ta[i].c.y*gpp[i].d.y + ta[i].c.z*gpp[i].e.y;
	sumoneindex.c.z = ta[i].c.x*gpp[i].c.z + ta[i].c.y*gpp[i].d.z + ta[i].c.z*gpp[i].e.z;
	sumoneindex.d.x = ta[i].d.x*gpp[i].c.x + ta[i].d.y*gpp[i].d.x + ta[i].d.z*gpp[i].e.x;
	sumoneindex.d.y = ta[i].d.x*gpp[i].c.y + ta[i].d.y*gpp[i].d.y + ta[i].d.z*gpp[i].e.y;
	sumoneindex.d.z = ta[i].d.x*gpp[i].c.z + ta[i].d.y*gpp[i].d.z + ta[i].d.z*gpp[i].e.z;
	sumtwoindex.x = sumoneindex.c.x*tb[i].c.x + sumoneindex.c.y*tb[i].c.y + sumoneindex.c.z*tb[i].c.z;
	sumtwoindex.y = sumoneindex.d.x*tb[i].d.x + sumoneindex.d.y*tb[i].d.y + sumoneindex.d.z*tb[i].d.z;
	float ts = -0.5f*(sumtwoindex.x+sumtwoindex.y);
	output[i].x = ts * (p1[i].x);
	output[i].y = ts * (p1[i].y);
}

__kernel void kernelradcontract_reswave_11(__global float44 * ta, __global float2 * p1, __global float44 * tb, __global float2 * q1, __global float44 * gpp, __global out float2 * output)
{
	uint i = get_global_id(0);
	float44 sumoneindex;
	float4 sumtwoindex;
	sumoneindex.c.x = ta[i].c.x*gpp[i].c.x + ta[i].c.y*gpp[i].d.x + ta[i].c.z*gpp[i].e.x;
	sumoneindex.c.y = ta[i].c.x*gpp[i].c.y + ta[i].c.y*gpp[i].d.y + ta[i].c.z*gpp[i].e.y;
	sumoneindex.c.z = ta[i].c.x*gpp[i].c.z + ta[i].c.y*gpp[i].d.z + ta[i].c.z*gpp[i].e.z;
	sumoneindex.d.x = ta[i].d.x*gpp[i].c.x + ta[i].d.y*gpp[i].d.x + ta[i].d.z*gpp[i].e.x;
	sumoneindex.d.y = ta[i].d.x*gpp[i].c.y + ta[i].d.y*gpp[i].d.y + ta[i].d.z*gpp[i].e.y;
	sumoneindex.d.z = ta[i].d.x*gpp[i].c.z + ta[i].d.y*gpp[i].d.z + ta[i].d.z*gpp[i].e.z;
	sumtwoindex.x = sumoneindex.c.x*tb[i].c.x + sumoneindex.c.y*tb[i].c.y + sumoneindex.c.z*tb[i].c.z;
	sumtwoindex.y = sumoneindex.d.x*tb[i].d.x + sumoneindex.d.y*tb[i].d.y + sumoneindex.d.z*tb[i].d.z;
	float ts = -0.5f*(sumtwoindex.x+sumtwoindex.y);
	output[i].x = ts * (p1[i].x*q1[i].x - p1[i].y*-q1[i].y);
	output[i].y = ts * (p1[i].x*-q1[i].y + p1[i].y*q1[i].x);
}

__kernel void kernelradcontract_reswave_12(__global float44 * ta, __global float2 * p1, __global float44 * tb, __global float2 * q1, __global float2 * q2, __global float44 * gpp, __global out float2 * output)
{
	uint i = get_global_id(0);
	float44 sumoneindex;
	float4 sumtwoindex;
	sumoneindex.c.x = ta[i].c.x*gpp[i].c.x + ta[i].c.y*gpp[i].d.x + ta[i].c.z*gpp[i].e.x;
	sumoneindex.c.y = ta[i].c.x*gpp[i].c.y + ta[i].c.y*gpp[i].d.y + ta[i].c.z*gpp[i].e.y;
	sumoneindex.c.z = ta[i].c.x*gpp[i].c.z + ta[i].c.y*gpp[i].d.z + ta[i].c.z*gpp[i].e.z;
	sumoneindex.d.x = ta[i].d.x*gpp[i].c.x + ta[i].d.y*gpp[i].d.x + ta[i].d.z*gpp[i].e.x;
	sumoneindex.d.y = ta[i].d.x*gpp[i].c.y + ta[i].d.y*gpp[i].d.y + ta[i].d.z*gpp[i].e.y;
	sumoneindex.d.z = ta[i].d.x*gpp[i].c.z + ta[i].d.y*gpp[i].d.z + ta[i].d.z*gpp[i].e.z;
	sumtwoindex.x = sumoneindex.c.x*tb[i].c.x + sumoneindex.c.y*tb[i].c.y + sumoneindex.c.z*tb[i].c.z;
	sumtwoindex.y = sumoneindex.d.x*tb[i].d.x + sumoneindex.d.y*tb[i].d.y + sumoneindex.d.z*tb[i].d.z;
	float ts = -0.5f*(sumtwoindex.x+sumtwoindex.y);
	output[i].x = ts * (p1[i].x*(q1[i].x*q2[i].x - q1[i].y*q2[i].y) - p1[i].y*-(q1[i].x*q2[i].y + q1[i].y*q2[i].x));
	output[i].y = ts * (p1[i].x*-(q1[i].x*q2[i].y + q1[i].y*q2[i].x) + p1[i].y*(q1[i].x*q2[i].x - q1[i].y*q2[i].y));
}

__kernel void kernelradcontract_reswave_20(__global float44 * ta, __global float2 * p1, __global float2 * p2, __global float44 * tb, __global float44 * gpp, __global out float2 * output)
{
	uint i = get_global_id(0);
	float44 sumoneindex;
	float4 sumtwoindex;
	sumoneindex.c.x = ta[i].c.x*gpp[i].c.x + ta[i].c.y*gpp[i].d.x + ta[i].c.z*gpp[i].e.x;
	sumoneindex.c.y = ta[i].c.x*gpp[i].c.y + ta[i].c.y*gpp[i].d.y + ta[i].c.z*gpp[i].e.y;
	sumoneindex.c.z = ta[i].c.x*gpp[i].c.z + ta[i].c.y*gpp[i].d.z + ta[i].c.z*gpp[i].e.z;
	sumoneindex.d.x = ta[i].d.x*gpp[i].c.x + ta[i].d.y*gpp[i].d.x + ta[i].d.z*gpp[i].e.x;
	sumoneindex.d.y = ta[i].d.x*gpp[i].c.y + ta[i].d.y*gpp[i].d.y + ta[i].d.z*gpp[i].e.y;
	sumoneindex.d.z = ta[i].d.x*gpp[i].c.z + ta[i].d.y*gpp[i].d.z + ta[i].d.z*gpp[i].e.z;
	sumtwoindex.x = sumoneindex.c.x*tb[i].c.x + sumoneindex.c.y*tb[i].c.y + sumoneindex.c.z*tb[i].c.z;
	sumtwoindex.y = sumoneindex.d.x*tb[i].d.x + sumoneindex.d.y*tb[i].d.y + sumoneindex.d.z*tb[i].d.z;
	float ts = -0.5f*(sumtwoindex.x+sumtwoindex.y);
	output[i].x = ts * (p1[i].x*p2[i].x - p1[i].y*p2[i].y);
	output[i].y = ts * (p1[i].x*p2[i].y + p1[i].y*p2[i].x);
}

__kernel void kernelradcontract_reswave_21(__global float44 * ta, __global float2 * p1, __global float2 * p2, __global float44 * tb, __global float2 * q1, __global float44 * gpp, __global out float2 * output)
{
	uint i = get_global_id(0);
	float44 sumoneindex;
	float4 sumtwoindex;
	sumoneindex.c.x = ta[i].c.x*gpp[i].c.x + ta[i].c.y*gpp[i].d.x + ta[i].c.z*gpp[i].e.x;
	sumoneindex.c.y = ta[i].c.x*gpp[i].c.y + ta[i].c.y*gpp[i].d.y + ta[i].c.z*gpp[i].e.y;
	sumoneindex.c.z = ta[i].c.x*gpp[i].c.z + ta[i].c.y*gpp[i].d.z + ta[i].c.z*gpp[i].e.z;
	sumoneindex.d.x = ta[i].d.x*gpp[i].c.x + ta[i].d.y*gpp[i].d.x + ta[i].d.z*gpp[i].e.x;
	sumoneindex.d.y = ta[i].d.x*gpp[i].c.y + ta[i].d.y*gpp[i].d.y + ta[i].d.z*gpp[i].e.y;
	sumoneindex.d.z = ta[i].d.x*gpp[i].c.z + ta[i].d.y*gpp[i].d.z + ta[i].d.z*gpp[i].e.z;
	sumtwoindex.x = sumoneindex.c.x*tb[i].c.x + sumoneindex.c.y*tb[i].c.y + sumoneindex.c.z*tb[i].c.z;
	sumtwoindex.y = sumoneindex.d.x*tb[i].d.x + sumoneindex.d.y*tb[i].d.y + sumoneindex.d.z*tb[i].d.z;
	float ts = -0.5f*(sumtwoindex.x+sumtwoindex.y);
	output[i].x = ts * ((p1[i].x*p2[i].x - p1[i].y*p2[i].y)*q1[i].x - (p1[i].x*p2[i].y + p1[i].y*p2[i].x)*-q1[i].y);
	output[i].y = ts * ((p1[i].x*p2[i].x - p1[i].y*p2[i].y)*-q1[i].y + (p1[i].x*p2[i].y + p1[i].y*p2[i].x)*q1[i].x);
}

__kernel void kernelradcontract_reswave_22(__global float44 * ta, __global float2 * p1, __global float2 * p2, __global float44 * tb, __global float2 * q1, __global float2 * q2, __global float44 * gpp, __global out float2 * output)
{
	uint i = get_global_id(0);
	float44 sumoneindex;
	float4 sumtwoindex;
	sumoneindex.c.x = ta[i].c.x*gpp[i].c.x + ta[i].c.y*gpp[i].d.x + ta[i].c.z*gpp[i].e.x;
	sumoneindex.c.y = ta[i].c.x*gpp[i].c.y + ta[i].c.y*gpp[i].d.y + ta[i].c.z*gpp[i].e.y;
	sumoneindex.c.z = ta[i].c.x*gpp[i].c.z + ta[i].c.y*gpp[i].d.z + ta[i].c.z*gpp[i].e.z;
	sumoneindex.d.x = ta[i].d.x*gpp[i].c.x + ta[i].d.y*gpp[i].d.x + ta[i].d.z*gpp[i].e.x;
	sumoneindex.d.y = ta[i].d.x*gpp[i].c.y + ta[i].d.y*gpp[i].d.y + ta[i].d.z*gpp[i].e.y;
	sumoneindex.d.z = ta[i].d.x*gpp[i].c.z + ta[i].d.y*gpp[i].d.z + ta[i].d.z*gpp[i].e.z;
	sumtwoindex.x = sumoneindex.c.x*tb[i].c.x + sumoneindex.c.y*tb[i].c.y + sumoneindex.c.z*tb[i].c.z;
	sumtwoindex.y = sumoneindex.d.x*tb[i].d.x + sumoneindex.d.y*tb[i].d.y + sumoneindex.d.z*tb[i].d.z;
	float ts = -0.5f*(sumtwoindex.x+sumtwoindex.y);
	output[i].x = ts * ((p1[i].x*p2[i].x - p1[i].y*p2[i].y)*(q1[i].x*q2[i].x - q1[i].y*q2[i].y) - (p1[i].x*p2[i].y + p1[i].y*p2[i].x)*-(q1[i].x*q2[i].y + q1[i].y*q2[i].x));
	output[i].y = ts * ((p1[i].x*p2[i].x - p1[i].y*p2[i].y)*-(q1[i].x*q2[i].y + q1[i].y*q2[i].x) + (p1[i].x*p2[i].y + p1[i].y*p2[i].x)*(q1[i].x*q2[i].x - q1[i].y*q2[i].y));
}
