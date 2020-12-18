#pragma once

#include <vector>
#include <TH1F.h>
#include <TH2F.h>
#include "TGraph.h"

using std::vector;
typedef vector<TH1F *> VTH1Fs;
typedef vector<VTH1Fs> VVTH1Fs;
typedef vector<TH2F *> VTH2Fs;
typedef vector<VTH2Fs> VVTH2Fs;
typedef vector<TGraph *> VTGraphs;

class GPUPlotset
{
public:
	GPUPlotset(void);
	/*##############################################################################
	# Function changed by LZH
	##############################################################################*/
	virtual ~GPUPlotset(void);
	void AddPlots(VTH1Fs plots) { mHistos.push_back(plots); };
	void AddPlots(VTH2Fs plots) { m2DHistos.push_back(plots); };
	void AddGraph(TGraph *graph) { mGraphs.push_back(graph); };
	void Format();
	void WriteRootfile(char *filename);
	/*##############################################################################
  	# Other function
  	##############################################################################*/
	void WritePsfile(char *filename, int xplots = 1, int yplots = 1);
	void WriteFiles(char *dirname, bool doeps = false, bool dopng = true);

protected:
	VVTH1Fs mHistos;
	VVTH2Fs m2DHistos;
	VTGraphs mGraphs;
};
