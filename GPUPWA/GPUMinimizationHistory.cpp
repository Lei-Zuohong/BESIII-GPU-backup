#include "GPUMinimizationHistory.h"
#include "GPUPartialWaveAnalysis.h"
#include <iostream>
#include <iomanip>

using std::left;
using std::right;
using std::setw;
string block = "****************************************************************************************************";

GPUMinimizationHistory::GPUMinimizationHistory(GPUPartialWaveAnalysis *ana,
											   const ROOT::Minuit2::MnUserParameters _InputParameters,
											   const ROOT::Minuit2::MnUserParameters _ResultParameters,
											   const double minimumLikelihood,
											   const double EDM,
											   const bool OK,
											   const int iterations) : mana(ana),
																	   m_InputParameters(_InputParameters),
																	   m_ResultParameters(_ResultParameters),
																	   m_minimumLikelihood(minimumLikelihood),
																	   mEDM(EDM),
																	   mOK(OK),
																	   mNiterations(iterations)
{
	mWaveNames = mana->GetActiveWaveNames();
	mMagnitudeParameters = mana->GetActiveMagnitudeParameters();
	mPhaseParameters = mana->GetActivePhaseParameters();
	mDynamicParameters = mana->GetActiveDynamicParameters();
	mNDynamicParameters = mana->GetNActiveDynamicParameters();
}

GPUMinimizationHistory::~GPUMinimizationHistory()
{
}

void GPUMinimizationHistory::Print(ostream &out)
{
	cout << block << endl;
	//PrintAnalysis(out);
	PrintFitInfo(out);
	//PrintWaves(out);
	//PrintInput(out);
	//PrintResult(out);
	//PrintInputRes(out);
	//PrintOutputRes(out);
	//PrintInputPara(out);
	//PrintOutputPara(out);
}

void GPUMinimizationHistory::PrintParameters(ROOT::Minuit2::MnUserParameters pars, ostream &out) const
{
	vector<double> pvec = pars.Params();
	out << "|" << setw(5) << left << "Index";
	out << "|" << setw(25) << left << "Parameter";
	out << "|" << setw(15) << left << "Value";
	out << "|" << setw(15) << left << "Error";
	out << "|" << endl;
	for (int i = 0; i < (int)pvec.size(); i++)
	{
		out << "|" << setw(5) << left << i;
		out << "|" << setw(25) << left << pars.Name(i);
		out << "|" << setw(15) << left << pvec[i];
		out << "|" << setw(15) << left << pars.Error(i);
		out << "|" << endl;
	}
}

void GPUMinimizationHistory::PrintInput(ostream &out) const
{
	out << setw(50) << right << "Parameters Input: " << endl;
	PrintParameters(m_InputParameters, out);
	out << block << endl;
}

void GPUMinimizationHistory::PrintResult(ostream &out) const
{
	out << setw(50) << right << "Parameters Output: " << endl;
	PrintParameters(m_ResultParameters, out);
	out << block << endl;
}

void GPUMinimizationHistory::PrintFitInfo(ostream &out) const
{
	out << setw(50) << right << "Status Report" << endl;
	out << "|" << setw(25) << left << "Fit";
	if (mOK)
		out << "|" << setw(15) << left << "Coverged";
	else
		out << "|" << setw(15) << left << "Failed";
	out << "|" << endl;
	out << "|" << setw(25) << left << "Minimum Likelihood";
	out << "|" << setw(15) << left << m_minimumLikelihood;
	out << "|" << endl;
	out << "|" << setw(25) << left << "Distance to Minimum";
	out << "|" << setw(15) << left << mEDM;
	out << "|" << endl;
	if (mNiterations)
	{
		out << "|" << setw(25) << left << "Iterations";
		out << "|" << setw(15) << left << mNiterations;
		out << "|" << endl;
	}
	cout << block << endl;
}

void GPUMinimizationHistory::PrintAnalysis(ostream &out)
{
	out << setw(50) << right << mana->GetName() << endl;
	out << "|" << setw(25) << left << "Data Events";
	out << "|" << setw(15) << left << mana->GetNumberData();
	out << "|" << endl;
	out << "|" << setw(25) << left << "MC Events Generated";
	out << "|" << setw(15) << left << mana->GetNumberMCGen();
	out << "|" << endl;
	out << "|" << setw(25) << left << "MC Events Accepted";
	out << "|" << setw(15) << left << mana->GetNumberMCAcc();
	out << "|" << endl;
	cout << block << endl;
}

void GPUMinimizationHistory::PrintWaves(ostream &out) const
{
	out << setw(50) << right << "Using waves:" << endl;
	out << "-------------------------------------------------------------------------------------------" << endl;
	out << "|" << setw(25) << left << "Wave Name";
	out << "|" << setw(15) << left << "Magn. in";
	out << "|" << setw(15) << left << "Magn. out";
	out << "|" << setw(15) << left << "Phase in";
	out << "|" << setw(15) << left << "Phase out";
	out << "|" << endl;
	out << "|" << setw(25) << left << "";
	out << "|" << setw(15) << left << "Dynamic in";
	out << "|" << setw(15) << left << "Dynamic out";
	out << "|" << setw(15) << left << "";
	out << "|" << setw(15) << left << "";
	out << "|" << endl;
	out << "-------------------------------------------------------------------------------------------" << endl;
	unsigned int dynparindex = 0;
	for (int i = 0; i < (int)mWaveNames.size(); i++)
	{
		out << "|" << setw(25) << left << mWaveNames[i];

		out << "|" << setw(15) << left << m_InputParameters.Value(mMagnitudeParameters[i]);
		if (m_InputParameters.Parameter(mMagnitudeParameters[i]).IsFixed())
			out << "|" << setw(15) << left << "(fixed)";
		else
			out << "|" << setw(15) << left << m_ResultParameters.Value(mMagnitudeParameters[i]);

		out << "|" << setw(15) << left << m_InputParameters.Value(mPhaseParameters[i]);
		if (m_InputParameters.Parameter(mPhaseParameters[i]).IsFixed())
			out << "|" << setw(15) << left << "(fixed)";
		else
			out << "|" << setw(15) << left << m_ResultParameters.Value(mPhaseParameters[i]);

		for (unsigned int n = 0; n < mNDynamicParameters[i]; n++)
		{
			out << "|" << endl;
			out << "|" << setw(25) << left << "";
			out << "|" << setw(15) << left << m_InputParameters.Value(*(mDynamicParameters[dynparindex]));
			if (m_InputParameters.Parameter(*(mDynamicParameters[dynparindex])).IsFixed())
				out << "|" << setw(15) << left << "(fixed)";
			else
				out << "|" << setw(15) << left << m_ResultParameters.Value(*(mDynamicParameters[dynparindex]));
			out << "|" << setw(15) << left << "";
			out << "|" << setw(15) << left << "";
			dynparindex++;
		}
		out << "|" << endl;
		out << "-------------------------------------------------------------------------------------------" << endl;
	}
	cout << block << endl;
}

void GPUMinimizationHistory::PrintParaParameters(ROOT::Minuit2::MnUserParameters pars, ostream &out) const
{
	int line1 = 25;
	int line2 = 15;
	int line3 = 15;
	int line4 = 15;
	int line5 = 15;
	vector<double> pvec = pars.Params();
	out << setw(line1) << left << "Parameter";
	out << "   ";
	out << setw(line2) << left << "Value";
	out << setw(line3) << left << "Error";
	out << setw(line4) << left << "Lower Limit";
	out << setw(line5) << left << "Upper Limit";
	out << endl;
	out << setw(line1) << left << pars.Name(0);
	out << " = ";
	out << setw(line2) << left << pvec[0];
	out << setw(line3) << left << pars.Error(0);
	out << setw(line4) << left << "0";
	out << setw(line5) << left << "500";
	out << endl;
	for (int i = 0; i < (int)mWaveNames.size(); i++)
	{
		out << setw(line1) << left << pars.Name(mMagnitudeParameters[i]);
		out << " = ";
		out << setw(line2) << left << pvec[mMagnitudeParameters[i]];
		out << setw(line3) << left << pars.Error(mMagnitudeParameters[i]);
		out << setw(line4) << left << (pars.Parameter(mMagnitudeParameters[i])).LowerLimit();
		out << setw(line5) << left << (pars.Parameter(mMagnitudeParameters[i])).UpperLimit();
		out << endl;
		out << setw(line1) << left << pars.Name(mPhaseParameters[i]);
		out << " = ";
		out << setw(line2) << left << pvec[mPhaseParameters[i]];
		out << setw(line3) << left << pars.Error(mPhaseParameters[i]);
		out << setw(line4) << left << (pars.Parameter(mPhaseParameters[i])).LowerLimit();
		out << setw(line5) << left << (pars.Parameter(mPhaseParameters[i])).UpperLimit();
		out << endl;
	}
}

void GPUMinimizationHistory::PrintResParameters(ROOT::Minuit2::MnUserParameters pars, ostream &out) const
{
	int line1 = 25;
	int line2 = 15;
	int line3 = 15;
	int line4 = 15;
	int line5 = 15;
	vector<double> pvec = pars.Params();
	out << setw(line1) << left << "Parameter";
	out << "   ";
	out << setw(line2) << left << "Value";
	out << setw(line3) << left << "Error";
	out << setw(line4) << left << "Lower Limit";
	out << setw(line5) << left << "Upper Limit";
	out << endl;
	for (int i = 1; i < (int)pvec.size(); i++)
	{
		if (fmod(i, 4) == 3)
		{
			out << setw(line1) << left << pars.Name(i);
			out << " = ";
			out << setw(line2) << left << pvec[i];
			out << setw(line3) << left << pars.Error(i);
			out << setw(line4) << left << pars.Parameter(i).LowerLimit();
			out << setw(line5) << left << pars.Parameter(i).UpperLimit();
			out << endl;
			out << setw(line1) << left << pars.Name(i + 1);
			out << " = ";
			out << setw(line2) << left << pvec[i + 1];
			out << setw(line3) << left << pars.Error(i + 1);
			out << setw(line4) << left << pars.Parameter(i + 1).LowerLimit();
			out << setw(line5) << left << pars.Parameter(i + 1).UpperLimit();
			out << endl;
		}
	}
}

void GPUMinimizationHistory::PrintInputRes(ostream &out) const
{
	PrintResParameters(m_InputParameters, out);
	cout << block << endl;
}

void GPUMinimizationHistory::PrintOutputRes(ostream &out) const
{
	PrintResParameters(m_ResultParameters, out);
	cout << block << endl;
}

void GPUMinimizationHistory::PrintInputPara(ostream &out) const
{
	PrintParaParameters(m_InputParameters, out);
	cout << block << endl;
}

void GPUMinimizationHistory::PrintOutputPara(ostream &out) const
{
	PrintParaParameters(m_ResultParameters, out);
	cout << block << endl;
}
