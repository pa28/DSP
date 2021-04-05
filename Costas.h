/**
 * @file Costas.h
 * @author Richard Buckley <richard.buckley@ieee.org>
 * @version 1.0
 * @date 2021-04-04
 */

#pragma once

#include <tuple>
#include "Iir.h"

namespace rose {

    /**
     * @class Costas
     * @brief
     */
    template<typename T>
    class Costas {
    protected:
        static constexpr std::size_t LowPassFilterOrder = 8;
        static constexpr std::size_t LoopFilterOrder = 8;

        using DataType = T;
        using LowPassFilterType = Iir::Butterworth::LowPass<LowPassFilterOrder>;
        using LoopFilterType = Iir::Butterworth::LowPass<LoopFilterOrder>;

        LowPassFilterType iLowPass{};
        LowPassFilterType qLowPass{};

        LoopFilterType loopFilter{};

        T sampleRate;
        T freqVCO;
        T freqLowPassCutOff;
        T dBLowPass;
        T freqLoopCutOff;
        T dbLoop;
        T dt{};
        T wVCO{};

//        T dt{}, wVCO{}, vcoI{}, vcoQ{}, bI{}, bQ{}, loopCtl{}, loopCorrection{};

    public:

        Costas() = delete;

        Costas(T sample_rate, T vco_freq, T low_pass_cut, T low_pass_db, T loop_cut, T loop_db) {
            sampleRate = sample_rate;
            freqVCO = vco_freq;
            freqLowPassCutOff = low_pass_cut;
            dBLowPass = low_pass_db;
            freqLoopCutOff = loop_cut;
            dbLoop = loop_db;
            dt = freqVCO / sampleRate;

            iLowPass.setup(sampleRate, freqLowPassCutOff);
            qLowPass.setup(sampleRate, freqLowPassCutOff);
            loopFilter.setup(sampleRate, freqLoopCutOff);
        }

        void sample(T a) {
            auto vcoI = sin(wVCO);
            auto vcoQ = cos(wVCO);

            auto aI = a * vcoI;
            auto aQ = a * vcoQ;

            auto bI = iLowPass.filter(vcoI);
            auto bQ = qLowPass.filter(vcoQ);

//            loopCtl = std::atan2(bI,bQ);

            printf( "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", a, vcoI, vcoQ, aI, aQ, bI, bQ, wVCO);
            wVCO = wVCO + (2. * M_PI * (dt));

            wVCO = fmod(wVCO, 2. * M_PI);
        }

//        std::tuple<T,T,T,T,T,T> state() const noexcept {
//            return std::make_tuple(wVCO, vcoI, vcoQ, bI, bQ, loopCtl);
//        }
    };
}

