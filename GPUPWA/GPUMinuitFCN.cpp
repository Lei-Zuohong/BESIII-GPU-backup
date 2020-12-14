#include "GPUMinuitFCN.h"

GPUMinuitFCN::GPUMinuitFCN(GPUPWACalculator *_calc) : mCalc(_calc)
{
}

GPUMinuitFCN::~GPUMinuitFCN(void)
{
}

vector<double> GPUMinuitFCN::Gradient(const vector<double> &x) const
{
	vector<double> a;
	mCalc->LikelihoodGradient(x, a);
	return a;
}

double GPUMinuitFCN::operator()(const vector<double> &x) const
{
	return mCalc->Likelihood(x);
}
