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
typedef std::vector<float> Vfloat;
typedef std::vector<double> Vdouble;

int cut_mc()
{
    // read
    TFile *tfile = new TFile("mc.root");
    TTree *ttree = tfile->Get("data");
    float E1, px1, py1, pz1;
    float E2, px2, py2, pz2;
    float E3, px3, py3, pz3;
    Vfloat VE1, Vpx1, Vpy1, Vpz1;
    Vfloat VE2, Vpx2, Vpy2, Vpz2;
    Vfloat VE3, Vpx3, Vpy3, Vpz3;
    ttree->SetBranchAddress("E1", &E1);
    ttree->SetBranchAddress("E2", &E2);
    ttree->SetBranchAddress("E3", &E3);
    ttree->SetBranchAddress("px1", &px1);
    ttree->SetBranchAddress("px2", &px2);
    ttree->SetBranchAddress("px3", &px3);
    ttree->SetBranchAddress("py1", &py1);
    ttree->SetBranchAddress("py2", &py2);
    ttree->SetBranchAddress("py3", &py3);
    ttree->SetBranchAddress("pz1", &pz1);
    ttree->SetBranchAddress("pz2", &pz2);
    ttree->SetBranchAddress("pz3", &pz3);
    for (int i = 0; i < 4; i++)
    {
        ttree->GetEntry(i);
        VE1.push_back(E1);
        VE2.push_back(E2);
        VE3.push_back(E3);
        Vpx1.push_back(px1);
        Vpx2.push_back(px2);
        Vpx3.push_back(px3);
        Vpy1.push_back(py1);
        Vpy2.push_back(py2);
        Vpy3.push_back(py3);
        Vpz1.push_back(pz1);
        Vpz2.push_back(pz2);
        Vpz3.push_back(pz3);
    }
    delete ttree;
    delete tfile;
    // write
    float nE1, npx1, npy1, npz1;
    float nE2, npx2, npy2, npz2;
    float nE3, npx3, npy3, npz3;
    TFile *tfile = new TFile("mc_4.root", "recreate");
    TTree *ttree = new TTree("data","data");
    TBranch *BE1 = ttree->Branch("E1", &nE1, "E1/F");
    TBranch *BE2 = ttree->Branch("E2", &nE2, "E2/F");
    TBranch *BE3 = ttree->Branch("E3", &nE3, "E3/F");
    TBranch *Bpx1 = ttree->Branch("px1", &npx1, "px1/F");
    TBranch *Bpx2 = ttree->Branch("px2", &npx2, "px2/F");
    TBranch *Bpx3 = ttree->Branch("px3", &npx3, "px3/F");
    TBranch *Bpy1 = ttree->Branch("py1", &npy1, "py1/F");
    TBranch *Bpy2 = ttree->Branch("py2", &npy2, "py2/F");
    TBranch *Bpy3 = ttree->Branch("py3", &npy3, "py3/F");
    TBranch *Bpz1 = ttree->Branch("pz1", &npz1, "pz1/F");
    TBranch *Bpz2 = ttree->Branch("pz2", &npz2, "pz2/F");
    TBranch *Bpz3 = ttree->Branch("pz3", &npz3, "pz3/F");
    for (int i = 0; i < 4; i++)
    {
        nE1 = VE1[i];
        nE2 = VE2[i];
        nE3 = VE3[i];
        npx1 = Vpx1[i];
        npx2 = Vpx2[i];
        npx3 = Vpx3[i];
        npy1 = Vpy1[i];
        npy2 = Vpy2[i];
        npy3 = Vpy3[i];
        npz1 = Vpz1[i];
        npz2 = Vpz2[i];
        npz3 = Vpz3[i];
        ttree->Fill();
    }
    ttree->Write();
    delete ttree;
    delete tfile;
    return 0;
}