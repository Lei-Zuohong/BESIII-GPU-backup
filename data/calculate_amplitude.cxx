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

#include "headc/bes_generator.h"

using namespace std;
typedef std::vector<float> Vfloat;

int calculate_amplitude()
{
    // read
    TFile *tfile = new TFile("data.root");
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
    for (int i = 0; i < 4; i++)
    {
        Vfloat v1 = bes_generator::new_Vfloat();
        v1[0] = Vpx1[i];
        v1[1] = Vpy1[i];
        v1[2] = Vpz1[i];
        v1[3] = VE1[i];
        Vfloat v2 = bes_generator::new_Vfloat();
        v2[0] = Vpx2[i];
        v2[1] = Vpy2[i];
        v2[2] = Vpz2[i];
        v2[3] = VE2[i];
        Vfloat v3 = bes_generator::new_Vfloat();
        v3[0] = Vpx3[i];
        v3[1] = Vpy3[i];
        v3[2] = Vpz3[i];
        v3[3] = VE3[i];
        Vfloat wave1 = bes_generator::new_Vfloat();
        Vfloat wave2 = bes_generator::new_Vfloat();
        Vfloat wave3 = bes_generator::new_Vfloat();
        Vfloat wave4 = bes_generator::new_Vfloat();
        Vfloat wave5 = bes_generator::new_Vfloat();
        wave1[0] = 1.000000;
        wave1[1] = 0.000000;
        wave2[0] = 0.800025;
        wave2[1] = -3.034010;
        wave3[0] = 0.308972;
        wave3[1] = 6.030600;
        wave4[0] = 0.004354;
        wave4[1] = 1.438700;
        wave5[0] = 0.102046;
        wave5[1] = 8.066140;
        float temp = bes_generator::amplitude(v1, v2, v3,
                                              wave1,
                                              wave2,
                                              wave3,
                                              wave4,
                                              wave5);
        cout << temp << endl;
    }
    return 0;
}