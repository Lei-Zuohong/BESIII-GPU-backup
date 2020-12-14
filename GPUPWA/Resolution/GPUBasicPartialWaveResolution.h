/// \file GPUBasicPartialWaveResolution.h

#pragma once
#ifndef GPUBASICPARTIALWAVERESOLUTION_H__
#define GPUBASICPARTIALWAVERESOLUTION_H__

#include "GPUPartialWave.h"
#include "GPUPropagatorType.h"
#include "GPUDataStream.h"
#include <cassert>
#include <vector>

class GPUBasicPartialWaveResolution: public GPUPartialWave
{
    public:
        GPUBasicPartialWaveResolution(char *_name, unsigned int nsets, GPUScalarPropagator &_propagator);
        virtual ~GPUBasicPartialWaveResolution(void);

        /// Get the list of the parameters, that control the dynamic part of this wave
        std::vector<unsigned int *> GetDynamicParameters() const { return mPropagator.GetParameters(); };
        /// Get the number the parameters, that control the dynamic part of this wave
        unsigned int GetNDynamicParameters() const { return mPropagator.GetNPars(); };

        /// Get the a parameter, that controls the dynamic part of this wave
        unsigned int GetDynamicParameter(unsigned int index) const { return mPropagator.GetParameter(index); };
        /// Set the number of the parameter, that controls the magnitude of this wave
        void SetMagnitudeParameter(int n) { SetParameter(0, n); };
        /// Set the number of the parameter, that controls the phase of this wave
        void SetPhaseParameter(int n) { SetParameter(1, n); };
        /// Set the numbers of the parameters, that control the dynamic part of this wave
        void SetDynamicParameters(std::vector<unsigned int *> n) { mPropagator.SetParameters(n); };
        /// Set the number of a parameter, that controls the dynamic part of this wave
        void SetDynamicParameter(unsigned int index, unsigned int parindex) {
            mPropagator.SetParameter(index, parindex);
        };

        /// Contract this wave with another one
        virtual GPUDataStream<float2> * Contract(GPUPartialWave * wave, unsigned int index, unsigned int block) = 0;

        /// Increase the Usecount for the caching mechanism
        virtual void IncreaseUsecount(unsigned int index) { mPropagator.IncreaseUsecount(index); };

    protected:
        char *mName;
        bool mActive;
        GPUScalarPropagator &mPropagator;
};

#endif
