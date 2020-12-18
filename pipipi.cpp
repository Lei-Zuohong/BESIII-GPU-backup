#pragma region 1.1 Include
#include "../GPUPWA/GPUStreamTensor.h"
#include "../GPUPWA/GPUComputedTensor.h"
#include "../GPUPWA/GPUComputedPropagator.h"
#include "../GPUPWA/GPUFactorizedPartialWave.h"
#include "../GPUPWA/GPUUnFactorizedPartialWave.h"
#include "../GPUPWA/GPUMetricTensor.h"
#include "../GPUPWA/GPUMetricTensorCMS.h"
#include "../GPUPWA/GPUAntisymmetricTensor.h"
#include "../GPUPWA/GPUOrbitalTensors.h"
#include "../GPUPWA/GPUPropagatorBreitWigner.h"
#include "../GPUPWA/GPUPropagatorFlatte2.h"
#include "../GPUPWA/GPUPropagatorSigma.h"
#include "../GPUPWA/GPUPropagatorMassDependentBreitWigner.h"
#include "../GPUPWA/GPUPartialWaveAnalysis.h"
#include "../GPUPWA/GPUPWAAmplitudeCalculator.h"
#include "../GPUPWA/GPUStreamInputVector.h"
#include "../GPUPWA/GPUStreamInputRootFileVector.h"
#include "../GPUPWA/GPUStreamInputTextFileVector.h"
#include "../GPUPWA/GPUPlotset.h"
#include "../GPUPWA/ConfigFile.h"
#include "../GPUPWA/ParaCfg.h"
#include "../GPUPWA/ResCfg.h"
#include "../GPUPWA/GPUChi2FitConstraint.h"
#include "../GPUPWA/GPUPropagatorGS.h"
#include "TFile.h"
#include "TRandom3.h"
#include "Minuit2/MnUserParameters.h"
#include <ctime>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <string>

#include "headc/bes_gpu.h"
#include "headc/bes_const.h"
#include "headc/common.h"
using namespace std;
#pragma endregion
#pragma region 1.2 constant
const double pi = bes_const::pi;
const double mpip = bes_const::particle_pipm.mass;
const double mpim = bes_const::particle_pipm.mass;
const double mpiz = bes_const::particle_piz.mass;
#pragma endregion
#define PLOT

int main(int argc, char *argv[])
{
#pragma region 2.0 Reading option
	common::OPTION_value myoption_value;
	common::OPTION_string myoption_string;
	myoption_value.readtxt("input_option_value.txt");
	myoption_string.readtxt("input_option_string.txt");
#pragma endregion
#pragma region 2.1 Reading file
	// 读取包文件
	GPUPartialWaveAnalysis *my_pwa = new GPUPartialWaveAnalysis((char *)"my_pwa", (char *)"input_files.txt", 2);
	my_pwa->SetMCIndex(1);
	GPUStreamInputRootFileVector &vector_pip = *new GPUStreamInputRootFileVector(my_pwa, my_pwa->GetDataFile(), "data", "px1", "py1", "pz1", "E1");
	GPUStreamInputRootFileVector &vector_pim = *new GPUStreamInputRootFileVector(my_pwa, my_pwa->GetDataFile(), "data", "px2", "py2", "pz2", "E2");
	GPUStreamInputRootFileVector &vector_piz = *new GPUStreamInputRootFileVector(my_pwa, my_pwa->GetDataFile(), "data", "px3", "py3", "pz3", "E3");
	vector_pip.SetFilename(my_pwa->GetMCFile(), my_pwa->GetMCIndex());
	vector_pim.SetFilename(my_pwa->GetMCFile(), my_pwa->GetMCIndex());
	vector_piz.SetFilename(my_pwa->GetMCFile(), my_pwa->GetMCIndex());
	// 读取root文件
	vector<int> index;
	vector<float> weight;
	index.push_back(int(myoption_value.get_value("number_data")));
	weight.push_back(1);
	index.push_back(int(myoption_value.get_value("number_sideband")));
	weight.push_back(-0.5);
	my_pwa->SetEventWeights(index, weight, 0);
	my_pwa->SetNumberMCGen(my_pwa->GetNumberMCAcc());
#pragma endregion
#pragma region 2.2 Construct unit
	GPUConstTensor4 &e_ijkl = *new GPUAntisymmetricTensor();

	GPUStreamVector &pip = vector_pip;
	GPUStreamVector &pim = vector_pim;
	GPUStreamVector &piz = vector_piz;
	GPUStreamVector &pipm = pip + pim;
	GPUStreamVector &pipz = pip + piz;
	GPUStreamVector &pimz = pim + piz;
	GPUStreamVector &pipipi = pip + pim + piz;

	GPUStreamScalar &pipm_m2 = pipm | pipm;
	GPUStreamScalar &pipz_m2 = pipz | pipz;
	GPUStreamScalar &pimz_m2 = pimz | pimz;
	GPUStreamScalar &pipm_m = sqrt(pipm_m2);
	GPUStreamScalar &pipz_m = sqrt(pipz_m2);
	GPUStreamScalar &pimz_m = sqrt(pimz_m2);

	GPUOrbitalTensors &q_all_rhop_pim = *new GPUOrbitalTensors(my_pwa, pipipi, pipz, pim);
	GPUOrbitalTensors &q_all_rhom_pip = *new GPUOrbitalTensors(my_pwa, pipipi, pimz, pip);
	GPUOrbitalTensors &q_all_rhoz_piz = *new GPUOrbitalTensors(my_pwa, pipipi, pipm, piz);
	GPUOrbitalTensors &q_rhop_pip_piz = *new GPUOrbitalTensors(my_pwa, pipz, pip, piz);
	GPUOrbitalTensors &q_rhom_pim_piz = *new GPUOrbitalTensors(my_pwa, pimz, pim, piz);
	GPUOrbitalTensors &q_rhoz_pip_pim = *new GPUOrbitalTensors(my_pwa, pipm, pip, pim);

	GPUStreamVector &t1_all_rhop_pim = q_all_rhop_pim.Spin1OrbitalTensor();
	GPUStreamVector &t1_all_rhom_pip = q_all_rhom_pip.Spin1OrbitalTensor();
	GPUStreamVector &t1_all_rhoz_piz = q_all_rhoz_piz.Spin1OrbitalTensor();
	GPUStreamVector &t1_rhop_pip_piz = q_rhop_pip_piz.Spin1OrbitalTensor();
	GPUStreamVector &t1_rhom_pim_piz = q_rhom_pim_piz.Spin1OrbitalTensor();
	GPUStreamVector &t1_rhoz_pip_pim = q_rhoz_pip_pim.Spin1OrbitalTensor();

	GPUStreamTensor3 &t3_all_rhop_pim = q_all_rhop_pim.Spin3OrbitalTensor();
	GPUStreamTensor3 &t3_all_rhom_pip = q_all_rhom_pip.Spin3OrbitalTensor();
	GPUStreamTensor3 &t3_all_rhoz_piz = q_all_rhoz_piz.Spin3OrbitalTensor();
	GPUStreamTensor3 &t3_rhop_pip_piz = q_rhop_pip_piz.Spin3OrbitalTensor();
	GPUStreamTensor3 &t3_rhom_pim_piz = q_rhom_pim_piz.Spin3OrbitalTensor();
	GPUStreamTensor3 &t3_rhoz_pip_pim = q_rhoz_pip_pim.Spin3OrbitalTensor();
#pragma endregion
#pragma region 2.3 Add wave
	// l = 1
	if (int(myoption_value.get_value("add_rho770pi")) == 1)
	{
		GPUPropagatorGS &propagator_rho770_p = *new GPUPropagatorGS((char *)"rho770", pipz_m2, mpip, mpiz, 3.0, 1);
		GPUPropagatorGS &propagator_rho770_m = *new GPUPropagatorGS((char *)"rho770", pimz_m2, mpim, mpiz, 3.0, 1);
		GPUPropagatorGS &propagator_rho770_z = *new GPUPropagatorGS((char *)"rho770", pipm_m2, mpip, mpim, 3.0, 1);
		GPUVectorPropagator &amp_rho770pi_p = (((e_ijkl | pipipi) | t1_rhop_pip_piz) | t1_all_rhop_pim) * propagator_rho770_p;
		GPUVectorPropagator &amp_rho770pi_m = (((e_ijkl | pipipi) | t1_rhom_pim_piz) | t1_all_rhom_pip) * propagator_rho770_m;
		GPUVectorPropagator &amp_rho770pi_z = (((e_ijkl | pipipi) | t1_rhoz_pip_pim) | t1_all_rhoz_piz) * propagator_rho770_z;
		GPUVectorPropagator &amp_rho770pi = amp_rho770pi_z - amp_rho770pi_p + amp_rho770pi_m;
		GPUUnFactorizedVectorPartialWave &wave_rho770pi = *new GPUUnFactorizedVectorPartialWave(amp_rho770pi, (char *)"wave_rho770pi");
		my_pwa->GetWaves()->AddPartialWave(wave_rho770pi);
	}
	if (int(myoption_value.get_value("add_rho1450pi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_rho1450_p = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1450", pipz_m2, 1, mpip, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rho1450_m = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1450", pimz_m2, 1, mpim, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rho1450_z = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1450", pipm_m2, 1, mpip, mpim);
		GPUVectorPropagator &amp_rho1450pi_p = (((e_ijkl | pipipi) | t1_rhop_pip_piz) | t1_all_rhop_pim) * propagator_rho1450_p;
		GPUVectorPropagator &amp_rho1450pi_m = (((e_ijkl | pipipi) | t1_rhom_pim_piz) | t1_all_rhom_pip) * propagator_rho1450_m;
		GPUVectorPropagator &amp_rho1450pi_z = (((e_ijkl | pipipi) | t1_rhoz_pip_pim) | t1_all_rhoz_piz) * propagator_rho1450_z;
		GPUVectorPropagator &amp_rho1450pi = amp_rho1450pi_z - amp_rho1450pi_p + amp_rho1450pi_m;
		GPUUnFactorizedVectorPartialWave &wave_rho1450pi = *new GPUUnFactorizedVectorPartialWave(amp_rho1450pi, (char *)"wave_rho1450pi");
		my_pwa->GetWaves()->AddPartialWave(wave_rho1450pi);
	}
	if (int(myoption_value.get_value("add_rho1570pi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_rho1570_p = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1570", pipz_m2, 1, mpip, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rho1570_m = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1570", pimz_m2, 1, mpim, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rho1570_z = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1570", pipm_m2, 1, mpip, mpim);
		GPUVectorPropagator &amp_rho1570pi_p = (((e_ijkl | pipipi) | t1_rhop_pip_piz) | t1_all_rhop_pim) * propagator_rho1570_p;
		GPUVectorPropagator &amp_rho1570pi_m = (((e_ijkl | pipipi) | t1_rhom_pim_piz) | t1_all_rhom_pip) * propagator_rho1570_m;
		GPUVectorPropagator &amp_rho1570pi_z = (((e_ijkl | pipipi) | t1_rhoz_pip_pim) | t1_all_rhoz_piz) * propagator_rho1570_z;
		GPUVectorPropagator &amp_rho1570pi = amp_rho1570pi_z - amp_rho1570pi_p + amp_rho1570pi_m;
		GPUUnFactorizedVectorPartialWave &wave_rho1570pi = *new GPUUnFactorizedVectorPartialWave(amp_rho1570pi, (char *)"wave_rho1570pi");
		my_pwa->GetWaves()->AddPartialWave(wave_rho1570pi);
	}
	if (int(myoption_value.get_value("add_rho1700pi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_rho1700_p = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1700", pipz_m2, 1, mpip, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rho1700_m = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1700", pimz_m2, 1, mpim, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rho1700_z = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1700", pipm_m2, 1, mpip, mpim);
		GPUVectorPropagator &amp_rho1700pi_p = (((e_ijkl | pipipi) | t1_rhop_pip_piz) | t1_all_rhop_pim) * propagator_rho1700_p;
		GPUVectorPropagator &amp_rho1700pi_m = (((e_ijkl | pipipi) | t1_rhom_pim_piz) | t1_all_rhom_pip) * propagator_rho1700_m;
		GPUVectorPropagator &amp_rho1700pi_z = (((e_ijkl | pipipi) | t1_rhoz_pip_pim) | t1_all_rhoz_piz) * propagator_rho1700_z;
		GPUVectorPropagator &amp_rho1700pi = amp_rho1700pi_z - amp_rho1700pi_p + amp_rho1700pi_m;
		GPUUnFactorizedVectorPartialWave &wave_rho1700pi = *new GPUUnFactorizedVectorPartialWave(amp_rho1700pi, (char *)"wave_rho1700pi");
		my_pwa->GetWaves()->AddPartialWave(wave_rho1700pi);
	}
	if (int(myoption_value.get_value("add_rhounknownpi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_rhounknown_p = *new GPUPropagatorMassDependentBreitWigner((char *)"rhounknown", pipz_m2, 1, mpip, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rhounknown_m = *new GPUPropagatorMassDependentBreitWigner((char *)"rhounknown", pimz_m2, 1, mpim, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rhounknown_z = *new GPUPropagatorMassDependentBreitWigner((char *)"rhounknown", pipm_m2, 1, mpip, mpim);
		GPUVectorPropagator &amp_rhounknownpi_p = (((e_ijkl | pipipi) | t1_rhop_pip_piz) | t1_all_rhop_pim) * propagator_rhounknown_p;
		GPUVectorPropagator &amp_rhounknownpi_m = (((e_ijkl | pipipi) | t1_rhom_pim_piz) | t1_all_rhom_pip) * propagator_rhounknown_m;
		GPUVectorPropagator &amp_rhounknownpi_z = (((e_ijkl | pipipi) | t1_rhoz_pip_pim) | t1_all_rhoz_piz) * propagator_rhounknown_z;
		GPUVectorPropagator &amp_rhounknownpi = amp_rhounknownpi_z - amp_rhounknownpi_p + amp_rhounknownpi_m;
		GPUUnFactorizedVectorPartialWave &wave_rhounknownpi = *new GPUUnFactorizedVectorPartialWave(amp_rhounknownpi, (char *)"wave_rhounknownpi");
		my_pwa->GetWaves()->AddPartialWave(wave_rhounknownpi);
	}
	if (int(myoption_value.get_value("add_omega782pi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_omega782 = *new GPUPropagatorMassDependentBreitWigner((char *)"omega782", pipm_m2, 1, mpip, mpim);
		GPUVectorPropagator &amp_omega782pi = (((e_ijkl | pipipi) | t1_rhoz_pip_pim) | t1_all_rhoz_piz) * propagator_omega782;
		GPUUnFactorizedVectorPartialWave &wave_omega782pi = *new GPUUnFactorizedVectorPartialWave(amp_omega782pi, (char *)"wave_omega782pi");
		my_pwa->GetWaves()->AddPartialWave(wave_omega782pi);
	}
	if (int(myoption_value.get_value("add_omega1420pi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_omega1420 = *new GPUPropagatorMassDependentBreitWigner((char *)"omega1420", pipm_m2, 1, mpip, mpim);
		GPUVectorPropagator &amp_omega1420pi = (((e_ijkl | pipipi) | t1_rhoz_pip_pim) | t1_all_rhoz_piz) * propagator_omega1420;
		GPUUnFactorizedVectorPartialWave &wave_omega1420pi = *new GPUUnFactorizedVectorPartialWave(amp_omega1420pi, (char *)"wave_omega1420pi");
		my_pwa->GetWaves()->AddPartialWave(wave_omega1420pi);
	}
	if (int(myoption_value.get_value("add_omega1650pi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_omega1650 = *new GPUPropagatorMassDependentBreitWigner((char *)"omega1650", pipm_m2, 1, mpip, mpim);
		GPUVectorPropagator &amp_omega1650pi = (((e_ijkl | pipipi) | t1_rhoz_pip_pim) | t1_all_rhoz_piz) * propagator_omega1650;
		GPUUnFactorizedVectorPartialWave &wave_omega1650pi = *new GPUUnFactorizedVectorPartialWave(amp_omega1650pi, (char *)"wave_omega1650pi");
		my_pwa->GetWaves()->AddPartialWave(wave_omega1650pi);
	}
	if (int(myoption_value.get_value("add_omegaunknownpi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_omegaunknown = *new GPUPropagatorMassDependentBreitWigner((char *)"omegaunknown", pipm_m2, 1, mpip, mpim);
		GPUVectorPropagator &amp_omegaunknownpi = (((e_ijkl | pipipi) | t1_rhoz_pip_pim) | t1_all_rhoz_piz) * propagator_omegaunknown;
		GPUUnFactorizedVectorPartialWave &wave_omegaunknownpi = *new GPUUnFactorizedVectorPartialWave(amp_omegaunknownpi, (char *)"wave_omegaunknownpi");
		my_pwa->GetWaves()->AddPartialWave(wave_omegaunknownpi);
	}
	// l = 3
	if (int(myoption_value.get_value("add_rho1690pi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_rho1690_p = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1690", pipz_m2, 3, mpip, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rho1690_m = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1690", pimz_m2, 3, mpim, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rho1690_z = *new GPUPropagatorMassDependentBreitWigner((char *)"rho1690", pipm_m2, 3, mpip, mpim);
		GPUVectorPropagator &amp_rho1690pi_p = (e_ijkl | pipipi) | (t3_all_rhop_pim || t3_rhop_pip_piz) * propagator_rho1690_p;
		GPUVectorPropagator &amp_rho1690pi_m = (e_ijkl | pipipi) | (t3_all_rhom_pip || t3_rhom_pim_piz) * propagator_rho1690_m;
		GPUVectorPropagator &amp_rho1690pi_z = (e_ijkl | pipipi) | (t3_all_rhoz_piz || t3_rhoz_pip_pim) * propagator_rho1690_z;
		GPUVectorPropagator &amp_rho1690pi = amp_rho1690pi_z - amp_rho1690pi_p + amp_rho1690pi_m;
		GPUUnFactorizedVectorPartialWave &wave_rho1690pi = *new GPUUnFactorizedVectorPartialWave(amp_rho1690pi, (char *)"wave_rho1690pi");
		my_pwa->GetWaves()->AddPartialWave(wave_rho1690pi);
	}
	if (int(myoption_value.get_value("add_rho3unknownpi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_rho3unknown_p = *new GPUPropagatorMassDependentBreitWigner((char *)"rho3unknown", pipz_m2, 3, mpip, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rho3unknown_m = *new GPUPropagatorMassDependentBreitWigner((char *)"rho3unknown", pimz_m2, 3, mpim, mpiz);
		GPUPropagatorMassDependentBreitWigner &propagator_rho3unknown_z = *new GPUPropagatorMassDependentBreitWigner((char *)"rho3unknown", pipm_m2, 3, mpip, mpim);
		GPUVectorPropagator &amp_rho3unknownpi_p = (e_ijkl | pipipi) | (t3_all_rhop_pim || t3_rhop_pip_piz) * propagator_rho3unknown_p;
		GPUVectorPropagator &amp_rho3unknownpi_m = (e_ijkl | pipipi) | (t3_all_rhom_pip || t3_rhom_pim_piz) * propagator_rho3unknown_m;
		GPUVectorPropagator &amp_rho3unknownpi_z = (e_ijkl | pipipi) | (t3_all_rhoz_piz || t3_rhoz_pip_pim) * propagator_rho3unknown_z;
		GPUVectorPropagator &amp_rho3unknownpi = amp_rho3unknownpi_z - amp_rho3unknownpi_p + amp_rho3unknownpi_m;
		GPUUnFactorizedVectorPartialWave &wave_rho3unknownpi = *new GPUUnFactorizedVectorPartialWave(amp_rho3unknownpi, (char *)"wave_rho3unknownpi");
		my_pwa->GetWaves()->AddPartialWave(wave_rho3unknownpi);
	}
	if (int(myoption_value.get_value("add_omega1670pi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_omega1670 = *new GPUPropagatorMassDependentBreitWigner((char *)"omega1670", pipm_m2, 3, mpip, mpim);
		GPUVectorPropagator &amp_omega1670pi = (e_ijkl | pipipi) | (t3_all_rhoz_piz || t3_rhoz_pip_pim) * propagator_omega1670;
		GPUUnFactorizedVectorPartialWave &wave_omega1670pi = *new GPUUnFactorizedVectorPartialWave(amp_omega1670pi, (char *)"wave_omega1670pi");
		my_pwa->GetWaves()->AddPartialWave(wave_omega1670pi);
	}
	if (int(myoption_value.get_value("add_omega3unknownpi")) == 1)
	{
		GPUPropagatorMassDependentBreitWigner &propagator_omega3unknown = *new GPUPropagatorMassDependentBreitWigner((char *)"omega3unknown", pipm_m2, 3, mpip, mpim);
		GPUVectorPropagator &amp_omega3unknownpi = (e_ijkl | pipipi) | (t3_all_rhoz_piz || t3_rhoz_pip_pim) * propagator_omega3unknown;
		GPUUnFactorizedVectorPartialWave &wave_omega3unknownpi = *new GPUUnFactorizedVectorPartialWave(amp_omega3unknownpi, (char *)"wave_omega3unknownpi");
		my_pwa->GetWaves()->AddPartialWave(wave_omega3unknownpi);
	}
#pragma endregion
#pragma region 3.1
	int number_waves = my_pwa->GetWaves()->GetNActiveWaves();
	vector_pip.ReadFile(0);
	vector_pim.ReadFile(0);
	vector_piz.ReadFile(0);
	vector_pip.ReadFile(1);
	vector_pim.ReadFile(1);
	vector_piz.ReadFile(1);
	my_pwa->MCIntegral();
	if (myoption_value.get_value("do_fit_minuit") == 1)
	{
		my_pwa->DoMultiFit(GPUPartialWaveAnalysis::MINUIT,
						   myoption_value.get_value("multifit_spread"),
						   myoption_value.get_value("multifit_times"),
						   false,
						   int(myoption_value.get_value("strategy_level")),
						   int(myoption_value.get_value("strategy_times")),
						   int(myoption_value.get_value("strategy_spread")));
	}
	my_pwa->InitCalculator();
	if (myoption_value.get_value("do_output_fraction") == 1)
	{
		double **ptx = my_pwa->GetPartialTotalXSection();
		char output_fraction[255];
		sprintf(output_fraction, "output_fraction.txt");
		ofstream output_fraction_file(output_fraction);
		bes_gpu::fraction_to_txt(ptx, number_waves, output_fraction_file);
		output_fraction_file.close();
	}
	if (myoption_value.get_value("do_output_amplitude") == 1)
	{
		my_pwa->GetAmplitudes((char *)"output_amplitude_data.txt",
							  (char *)"output_amplitude_mc.txt");
	}
#pragma endregion
#pragma region 3.2
	if (myoption_value.get_value("do_output_root") == 1)
	{
		my_pwa->Reset(0);
		float **dcs = my_pwa->GetMCDcs(true);
		float *pweights = my_pwa->GetEventWeights();
		my_pwa->Reset(1);
		GPUStreamScalar &pip_a = costheta(pip);
		GPUStreamScalar &pim_a = costheta(pim);
		GPUStreamScalar &piz_a = costheta(piz);
		GPUPlotset *plotset = new GPUPlotset();
		plotset->AddPlots(pipm_m.Plot((char *)"pipm_m", (char *)"M(#pi^{+}#pi^{-})(GeV/c^{2}); Events/10 MeV/c^{2}", 90, 0.2, 2., dcs, number_waves, true, pweights));
		plotset->AddPlots(pipz_m.Plot((char *)"pipz_m", (char *)"M(#pi^{+}#pi^{0})(GeV/c^{2}); Events/10 MeV/c^{2}", 90, 0.2, 2., dcs, number_waves, true, pweights));
		plotset->AddPlots(pimz_m.Plot((char *)"pimz_m", (char *)"M(#pi^{-}#pi^{0})(GeV/c^{2}); Events/10 MeV/c^{2}", 90, 0.2, 2., dcs, number_waves, true, pweights));
		plotset->AddPlots(pip_a.Plot((char *)"pip_a", (char *)"cos theta of #pi^{+} in ECMS; cos(#theta_{#pi^{+}}) ", 50, -1, 1, dcs, number_waves, true, pweights));
		plotset->AddPlots(pim_a.Plot((char *)"pim_a", (char *)"cos theta of #pi^{-} in ECMS; cos(#theta_{#pi^{-}}) ", 50, -1, 1, dcs, number_waves, true, pweights));
		plotset->AddPlots(piz_a.Plot((char *)"piz_a", (char *)"cos theta of #pi^{0} in ECMS; cos(#theta_{#pi^{0}}) ", 50, -1, 1, dcs, number_waves, true, pweights));
		plotset->Format();
		plotset->WriteRootfile((char *)(myoption_string.get_value("root_output").data()));
	}
#pragma endregion
	return 0;
}
