#pragma once

#include <iostream>
#include <vector>

#include "Minuit2/MnUserParameters.h"

#include "GPURuncounter.h"
#include "GPUMinimizationHistory.h"
#include "GPUDataStream.h"
#include "GPUDataDependentObjectList.h"
#include "GPUSetOfPartialWaves.h"
#include "GPUPartialWave.h"

#include "Status.h"
#include "ParaCfg.h"
#include "ResCfg.h"
#include "ConfigFile.h"

class GPUPWACalculator;
class GPUFitConstraintList;
class GPUFitConstraint;
class GPULookupTable;
class TGraph;

using std::cout;
using std::endl;
using std::vector;
typedef std::vector<double> Vdouble;
typedef std::vector<Vdouble> VVdouble;

class GPUPartialWaveAnalysis : public GPUDataDependentObjectList
{
public:
	enum FITTER
	{
		FUMILI,
		MINUIT,
		MINUITGRAD,
		MINUITMINOS,
		OLDFUMILI,
		INITONLY
	};

	/// Constructor
	/**
	\param _name: This name is used for recognising files belonging to this analysis and related administrative things.
	\param _filename: This is the name of a file that contains the names of the resonance and fit parameter files to be read
	\param _ndatasets: Number of data/MC sets used in this analysis - default is 2 for Data at index 0 and MC at index 1
	**/
	GPUPartialWaveAnalysis(char *_name, char *_filename, int _ndatasets = 2);
	/// Destructor
	virtual ~GPUPartialWaveAnalysis(void);
	/// Get the name of this analysis
	char *GetName() const { return mName; };
	/// Get the name of this analysis with spaces stripped
	char *GetNoSpaceName() const { return mNoSpaceName; };
	// Handling of Dataset parameters
	/// Set the index of the data streams (default is 0)
	void SetDataIndex(unsigned int _i = 0) { mDataIndex = _i; };
	/// Set the index of the MC streams (default is 1)
	void SetMCIndex(unsigned int _i = 1) { mMCIndex = _i; };
	/// Get the index of the data streams
	unsigned int GetDataIndex() const { return mDataIndex; };
	/// Get the index of the MC streams
	unsigned int GetMCIndex() const { return mMCIndex; };
	/// Get the set of partial waves
	GPUSetOfPartialWaves *GetWaves() const { return mWaves; };
	// Mapping of fit parameters to waves
	/// Set the magnitude parameter number of a wave
	void SetMagnitudeParameter(int wave, int parnum)
	{
		(*mWaves)[wave]->SetMagnitudeParameter(parnum);
	};
	/// Set the phase parameter number of a wave
	void SetPhaseParameter(int wave, int parnum) { (*mWaves)[wave]->SetPhaseParameter(parnum); };
	// Handling the fit parameters - interface to MnUserParameters provided for convenience
	/// Get a pointer to the fit parameters
	ROOT::Minuit2::MnUserParameters *GetParameters() const { return mParameters; };
	/// Get A copy of the fit parameters
	ROOT::Minuit2::MnUserParameters *GetParameterCopy() const
	{
		return new ROOT::Minuit2::MnUserParameters(*mParameters);
	};
	/// Set the fit parameters
	void SetParameters(ROOT::Minuit2::MnUserParameters *pars)
	{
		assert(pars);
		mParameters = pars;
	};
	/// Get the vector of fit parameters
	vector<double> GetParameterVector() const { return mParameters->Params(); };
	/// Get the i-th fit parameters
	double GetParameter(int i) const { return mParameters->Value(i); };
	/// Get the value of a fit parameter known by name
	double GetParameter(const char *s) const { return mParameters->Value(mParameters->Index(s)); };
	/// Get the index of a fit parameter known by name
	unsigned int GetParameterIndex(const char *s) const { return mParameters->Index(s); };
	/// Get the Name of a parameter
	const char *GetParameterName(int i) const { return mParameters->Name(i); };
	/// Set the vector of fit parameters
	void SetParameterVector(vector<double> _vec) const;
	/// Set the i-th fit parameter
	void SetParameter(int i, double val) { mParameters->SetValue(i, val); };
	/// Set the i-th fit parameter from a ParaCfg object (usually read from a ConfigFile)
	void SetParameter(int i, ParaCfg cfg)
	{
		SetParameter(i, cfg.v);
		if (cfg.e > 0)
		{
			SetParameterError(i, cfg.e);
		}
		else
		{
			FixParameter(i, cfg.v);
			SetParameterError(i, -cfg.e);
		}
		if (cfg.l != 999 && cfg.u != 999)
		{
			LimitParameter(i, cfg.l, cfg.u);
		}
		else if (cfg.l != 999)
		{
			LimitParameterLow(i, cfg.l);
		}
		else if (cfg.u != 999)
		{
			LimitParameterHigh(i, cfg.u);
		}
	}
	/// Fix a paramater
	void FixParameter(int i, double val)
	{
		mParameters->SetValue(i, val);
		mParameters->Fix(i);
	};
	/// Fix a parameter to its current value
	void FixParameter(int i)
	{
		mParameters->Fix(i);
	};
	/// Free a parameter
	void ReleaseParameter(int i) { mParameters->Release(i); };
	/// Limit a parameter (two-sided)
	void LimitParameter(int i, double low, double high) { mParameters->SetLimits(i, low, high); };
	/// Limit a parameter (lower limit)
	void LimitParameterLow(int i, double low) { mParameters->SetLowerLimit(i, low); };
	/// Limit a parameter (upper limit)
	void LimitParameterHigh(int i, double high) { mParameters->SetUpperLimit(i, high); };
	/// Add a parameter
	int AddParameter(char *name, double _val, double _err)
	{
		mParameters->Add(name, _val, _err);
		return mParameters->Index(name);
	};
	/// Set Parameter Error
	void SetParameterError(int i, double err) { mParameters->SetError(i, err); };
	/// Set Parameter ErrorVector
	void SetParameterErrorVector(std::vector<double> errs)
	{
		for (int i = 0; i < (int)errs.size(); i++)
			mParameters->SetError(i, errs[i]);
	};
	/// Print the parametres to stdout
	void PrintParameters() const;

	// Handling event weights
	/// Set event weights from an array of weights
	/** It is the responsibility of the user that the weights array is at least as long
	 as the number of data events; This function should only be called after the data
	files are set, as it needs to know the number of data events**/
	void SetEventWeights(float *weights, unsigned int index = 0);
	/// Set event weights to a constant
	/** This function should only be called after the data
	    files are set, as it needs to know the number of data events
	**/
	void SetEventWeights(float weight = 1.0, unsigned int index = 0);
	/// Set event weights from a vector of event numbers and a vector of weights
	/** The two vector should contain the same number of entries. Using this function,
	 the first nums[0] events will be assigned a weight of weights[0] and so on.
	If the sum of nums is different from the number of data events, this will in a program abort.
	.This function should only be called after the data
	files are set, as it needs to know the number of data events **/
	void SetEventWeights(vector<int> nums, vector<float> weights, unsigned int index = 0);

	/// Get the event weights (CPU)
	float *GetEventWeights(unsigned int index = 0) { return mCPUWeights[index]; };

	/// Get the event weights (GPU)
	GPUDataStream<float4> **GetWeights(unsigned int index = 0) { return mWeights[index]; };

	/// Get the event weights (GPU), floats
	GPUDataStream<float> **GetWeightsSingle(unsigned int index = 0) { return mWeightsSingle[index]; };

	// Handling event numbers
	/// Set Number of generated MC events
	void SetNumberMCGen(int nmc) { mNMCgen = nmc; };
	/// Get Number of generated MC events
	int GetNumberMCGen() const { return mNMCgen; };
	/// Get Number of accepted MC events
	/* Only returns a meaningful value after the MC input files have been set */
	int GetNumberMCAcc()
	{
		if (GetNevents(GetMCIndex()) == -1)
			cout << "WARNING: Set MC files before calling GetNumberMCAcc()!" << endl;
		return GetNevents(GetMCIndex());
	};
	/// Get Number of accepted MC events
	/* Only returns a meaningful value after the MC input files have been set */
	double GetSumOfWeights(int index)
	{
		if (!msumofweights[index])
		{
			msumofweights[index] = GetNevents(index);
		}
		return msumofweights[index];
	};
	/// Get Number of Data events
	/* Only returns a meaningful value after the data input files have been set */
	int GetNumberData()
	{
		if (GetNevents(GetDataIndex()) == -1)
		{
			cout << "WARNING: Set Data files before calling GetNumberData()!" << endl;
		}
		return GetNevents(GetDataIndex());
	};

	// Handle constraints
	/// Add a constraint to this partial wave analysis
	void AddConstraint(GPUFitConstraint *constraint);
	/// Remove a constraint from this partial wave analysis
	void RemoveConstraint(unsigned int index);
	/// Get a pointer to a constraint
	GPUFitConstraint *GetConstraint(unsigned int index) const;
	/// Get the constraint list
	GPUFitConstraintList *GetConstraintList() const { return mConstraintList; };
	/// Get the log likelihood contribution of the constraints
	double GetConstraintLHContribution() const;
	/// Get the log likelihood gradient contribution of the constraints with regard to the parameter at index
	double GetConstraintGradientContribution(unsigned int index) const;
	/// Get the log likelihood hessian contribution (in the FUmili approximation) of the constraints with regard to the parameters i,j
	double GetConstraintHessianContribution(unsigned int i, unsigned int j) const;
	/// Describe constraints
	void DescribeConstraints(std::ostream &outstream) const;
	/// Report constraints contribution to likelihood
	void ReportConstraints(std::ostream &outstream) const;
	/// Construct a lookup table for an index
	GPULookupTable *GetLookupTable(unsigned int index);
	// The serious stuff
	/// Calculate the MC integral coefficients for this partial wave analysis
	/** /param all: Select wheter the parameters should be calculated for all (recommended) or only the active waves
	    On calling this function, all input streams should be linked to the correct MC data files
	    If it is needed in the contraction (radiative decay to mesons), g_mu_nu^(perp perp) has to be set
	**/
	virtual void MCIntegral(bool all = true);
	/// Perform the PWA fit
	virtual STATUS DoFit(FITTER fitter);

	/// Write the amplitudes for the current parameter set to a file
	virtual void GetAmplitudes(char *filename_data, char *filename_mc);

	/// Perform a set of fits with different sets of input values
	virtual STATUS DoDynamicFit(FITTER fitter,
								double deltas,
								int ngood,
								VVdouble inputs,
								bool stopatconvergence);
	/// Perform a set of fits with different sets of input values
	virtual STATUS DoMultiFit(FITTER fitter,
							  VVdouble inputs,
							  bool stopatconvergence,
							  int strategy_level, int strategy_times, int strategy_spread);
	/// Perform a set of fits with different input values chosen at random
	/** All free parameters are varied simultaneously according to a flat pdf
	 *  within spread times the error specified for the respective parameter
	 *  in total nfit configurations are generated and tested. If the fit converges
	 *  to a likelihood value within deltas of the minimum likelihood value ngood
	 *  times, no further fits are performed. (Written by Xiaoling Li)
	 */
	virtual STATUS DoDynamicFit(FITTER fitter,					 ///< fitter to be used
								double deltas = 0.5,			 ///< Maximum difference to minimum for stop
								int ngood = 5,					 ///< Number of times fit has to converge within deltas
								double spread = 5.0,			 ///< spread (in parameter errors) for the variations
								int nfit = 1000,				 ///< Number of configurations to generate and fit
								bool stopatconvergence = false); ///< Stop after the first fit has converged
	virtual STATUS DoMultiFit(FITTER fitter,					 ///< fitter to be used
							  double spread = 5.0,				 ///< spread (in parameter errors) for the variations
							  int nfit = 10,					 ///< Number of configurations to generate and fit
							  bool stopatconvergence = false,	 ///< Stop after the first fit has converged
							  int strategy_level = 1,
							  int strategy_times = 100000,
							  int strategy_spread = 20);
	/// Perform a parameter scan
	virtual TGraph *ScanParameter(FITTER fitter,	  ///< fitter to be used
								  int parameterindex, ///< Parameter to be scanned
								  double min,		  ///< lowest value for the parameter
								  double max,		  ///< highest value for the parameter
								  int steps);		  ///< Number of steps to test
	/// Create the matrix encoding parameter-wave relations with fixed masses and withs
	virtual float **parmatrix() const;
	/// Create the matrix encoding parameter-wave relations with free masses and widths
	virtual float **allparmatrix() const;
	/// create a stream from the amplitude parameters used for the active waves, using cartesian complex numbers
	virtual GPUDataStream<float2> *pars2streamCartesian() const;
	/// create a stream from the amplitude parameters used for the active waves, using polar complex numbers
	virtual GPUDataStream<float2> *pars2streamPolar() const;
	/// create a stream from the amplitude and resonance parameters used for the active waves, using cartesian complex numbers
	virtual GPUDataStream<float4> *allpars2streamCartesian() const;
	/// create a stream from the amplitude and resonance parameters used for the active waves, using polar complex numbers
	virtual GPUDataStream<float4> *allpars2streamPolar() const;
	/// get minimum Likelihood value
	/** If not fit has been performed, the returned value does not correspond to any Likelihood **/
	virtual double GetMinLikelihood() const { return mLikelihood; };
	/// Get an array of weight values for the MC events in order to produce plots
	virtual float **GetMCDcs(bool dooffdiagonal = false, ///< Do or do not compute interference terms
							 int nblocks = -1			 ///< Maximum number ob blocks to plot for MC (-1 = all)
	);
	/// Get The calculator object
	/** The returned pointer will only be valid if you called DoFit(INITONLY) beforehand.
	 It can be used for debugging or calculating Likelihoods at specific points**/
	virtual GPUPWACalculator *GetCalculator() const { return mcalc; };
	/// Get parameter file
	virtual const ConfigFile *GetParameterFile() const { return mParaFile; };
	/// Get resonance file
	virtual const ConfigFile *GetResonanceFile() const { return mResFile; };
	/// Get steering file
	virtual const ConfigFile *GetSteeringFile() const { return mSteeringFile; };
	/// Get the data file
	virtual const char *GetDataFile() const { return mInputFiles[mDataIndex]; };
	/// Get the MC file
	virtual const char *GetMCFile() const { return mInputFiles[mMCIndex]; };
	/// Get a MC file
	virtual const char *GetMCFile(unsigned int index) const { return mInputFiles[index]; };
	/// Get the number of partial waves
	virtual int GetNPartialWaves() const { return mWaves->GetNWaves(); };
	/// Get the number of partial waves
	virtual int GetNActivePartialWaves() const { return mWaves->GetNActiveWaves(); };
	/// Get a partial wave
	virtual GPUPartialWave *GetPartialWave(int n) const { return mWaves->GetPartialWave(n); }
	/// Get a vector of wave names
	virtual std::vector<char *> GetWaveNames() const { return mWaves->GetWaveNames(); };
	/// Get a vector of wave names for the active waves
	virtual std::vector<char *> GetActiveWaveNames() const { return mWaves->GetWaveNames(); };
	/// Get a vector of the magnitude parameter indices
	virtual std::vector<unsigned int> GetMagnitudeParameters() const { return mWaves->GetMagnitudeParameters(); };
	/// Get a vector of the phase parameter indices
	virtual std::vector<unsigned int> GetPhaseParameters() const { return mWaves->GetPhaseParameters(); };
	/// Get a vector off the dynamic parameter indices
	virtual std::vector<unsigned int *> GetDynamicParameters() const { return mWaves->GetDynamicParameters(); };
	/// Get a vector off the number of dynamic parameter indices
	virtual std::vector<unsigned int> GetNDynamicParameters() const { return mWaves->GetNDynamicParameters(); };
	/// Get a vector of the magnitude parameter indices for active waves
	virtual std::vector<unsigned int> GetActiveMagnitudeParameters() const { return mWaves->GetActiveMagnitudeParameters(); };
	/// Get a vector of the phase parameter indices for active waves
	virtual std::vector<unsigned int> GetActivePhaseParameters() const { return mWaves->GetActivePhaseParameters(); };
	/// Get a vector of the dynamic parameter indices for active waves
	virtual std::vector<unsigned int *> GetActiveDynamicParameters() const { return mWaves->GetActiveDynamicParameters(); };
	/// Get a vector of number of the dynamic parameter indices for active waves
	virtual std::vector<unsigned int> GetNActiveDynamicParameters() const { return mWaves->GetNActiveDynamicParameters(); };
	/// Set the usecount (used for cache handling)
	virtual void SetUsecount(int index) { mWaves->SetUsecount(index); };
	/// Set the usecount for the orbital part only (used for cache handling)
	virtual void SetUsecountOrbital(int index) { mWaves->SetUsecountOrbital(index); };

	/// Check whether all resonance parameters are fixed or a free parameter fit has to be attempted
	virtual bool HasFreeResonanceParameters() const;

	/// Get the total cross section of each partial wave
	virtual double **GetPartialTotalXSection() const;
	/// Get the derivation of the total cross section of each partial wave
	virtual double **GetPartialTotalXSectionDerivative(int parnum) const;
	/// initalize GPUPWACalculator::mcalc, used by GPUChi2FrFitConstraint
	virtual void InitCalculator();

	/// Set lookuptables - they will be used instead of the partial wave set
	virtual void SetLookupTables(GPULookupTable *datatable, GPULookupTable *mctable)
	{
		mLookupData = datatable;
		mLookupMC = mctable;
	}

protected:
	// 拟合方法
	virtual STATUS MinuitMinimize(GPUPWACalculator *calc, GPUMinimizationHistory *&minhist,
								  int strategy_level = 1, int strategy_times = 100000, int strategy_spread = 20);
	virtual STATUS MinuitGradMinimize(GPUPWACalculator *calc, GPUMinimizationHistory *&minhist);
	virtual STATUS FumiliMinimize(GPUPWACalculator *calc, GPUMinimizationHistory *&minhist);
	virtual STATUS OldfumiliMinimize(GPUPWACalculator *calc, GPUMinimizationHistory *&minhist);
	virtual STATUS MinuitMinosMinimize(GPUPWACalculator *calc, GPUMinimizationHistory *&minhist);
	// 通用变量
	char *mName;
	char *mNoSpaceName;
	char *mFilename;
	std::vector<char *> mInputFiles;
	ConfigFile *mParaFile;
	ConfigFile *mResFile;
	ConfigFile *mSteeringFile;
	int mNMCgen;
	unsigned int mDataIndex;
	unsigned int mMCIndex;
	std::vector<GPUDataStream<float4> **> mWeights;
	std::vector<GPUDataStream<float> **> mWeightsSingle;
	std::vector<float *> mCPUWeights;
	std::vector<double> msumofweights;
	// 拟合变量
	mutable ROOT::Minuit2::MnUserParameters *mParameters;
	GPUPWACalculator *mcalc;
	double mLikelihood;
	GPURuncounter *mRuncounter;
	mutable GPUDataStream<float2> *mparscartesian;
	mutable GPUDataStream<float2> *mparspolar;
	GPUMinimizationHistory *mLastFitHistory;
	GPUMinimizationHistory *mBestFitHistory;
	GPUSetOfPartialWaves *mWaves;
	GPULookupTable *mLookupData;
	GPULookupTable *mLookupMC;
	GPUFitConstraintList *mConstraintList;
};
