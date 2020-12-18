#include <TCanvas.h>
#include <TH2.h>
#include "TTree.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <vector>
using namespace std;
typedef std::vector<double> Vdouble;

int test()
{
    cout << "this is test" << endl;
    Vdouble pipm_m;
    Vdouble weight;
    TFile *tfile = new TFile("data.root");
    TTree *ttree = tfile->Get("data");
    float E1, E2, px1, px2, py1, py2, pz1, pz2;
    ttree->SetBranchAddress("E1", &E1);
    ttree->SetBranchAddress("E2", &E2);
    ttree->SetBranchAddress("px1", &px1);
    ttree->SetBranchAddress("px2", &px2);
    ttree->SetBranchAddress("py1", &py1);
    ttree->SetBranchAddress("py2", &py2);
    ttree->SetBranchAddress("pz1", &pz1);
    ttree->SetBranchAddress("pz2", &pz2);
    float temp;
    for (int i = 0; i < 9724; i++)
    {
        ttree->GetEntry(i);
        temp = pow(E1 + E2, 2) - pow(px1 + px2, 2) - pow(py1 + py2, 2) - pow(py1 + py2, 2);
        pipm_m.push_back(sqrt(temp));
    }
    delete ttree;
    delete tfile;

    ifstream myfile("../output_amplitude_data.txt");
    string temp1;
    float temp2;
    for (int i = 0; i < 9724; i++)
    {
        myfile >> temp1;
        temp2 = atof(temp1.c_str());
        weight.push_back(temp2);
    }

    TFile *tfile = new TFile("hist.root", "recreate");
    TH1F *thist = new TH1F("hist", "hist", 250, 0.0, 2.5);
    for (int i = 0; i < 9724; i++)
    {
        thist->Fill(pipm_m[i], 1/weight[i]);
    }
    thist->Write();
    delete thist;
    delete tfile;

    return 0;
}