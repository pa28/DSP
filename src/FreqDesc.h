/**
 * @file FreqDesc.h
 * @author Richard Buckley <richard.buckley@ieee.org>
 * @version 1.0
 * @date 2021-04-05
 */

#pragma once

#include "Iir.h"

namespace rose {

    /**
     * @class FreqDesc
     * @brief A frequency discriminator.
     */
    struct FreqDesc {
        Iir::Butterworth::LowPass<2> lpfI;  ///< Low pass filter for post I mixer filtering.
        Iir::Butterworth::LowPass<2> lpfQ;  ///< Low pass filter for post Q mixer filtering.

        double dt,      ///< delta time parameter for VCO free-running frequency
        w,       ///< the current phase of the VCO.
        phi,     ///< The last phase output by the quadrature mixer.
        bPhi,    ///< The output of the loop control filter.
        aI,      ///< Input signal after mixing with I phase of the VCO.
        aQ,      ///< Input signal after mixing with Q phase of the VCO.
        bI,      ///< Filtered aI
        bQ;      ///< Filtered aQ

        double fSamp{}; ///< Sample rate to convert results to signal related values.

        FreqDesc() = delete;

        /**
         * @brief Construct and configure the frequency discriminator.
         * @param sampleRate The input signal sampling rate in Hz.
         * @param centreFrequency The local oscillator frequency in Hz.
         */
        FreqDesc(double sampleRate, double centreFrequency, double modulationRate) :
                lpfI(), lpfQ() {
            fSamp = sampleRate;
            lpfI.setup(sampleRate, modulationRate);
            lpfQ.setup(sampleRate, modulationRate);
            dt = centreFrequency / sampleRate;
            w = 0.;
            phi = bPhi = aI = aQ = bI = bQ = 0.;
        }

        /**
         * @brief Process an input sample.
         * @param a The input sample.
         * @return The frequency deviation from the local oscillator.
         */
        double sample(double a) {
            w = std::fmod(w + (2. * M_PI * (dt)), 2. * M_PI);
            aI = std::sin(w) * a;
            aQ = std::cos(w) * a;

            bI = lpfI.filter(aI);
            bQ = lpfQ.filter(aQ);

            phi = std::atan2(bQ, bI);

            auto r = phi - bPhi;
            if (r < - M_PI) r += 2. * M_PI;
            if (r > M_PI) r -= 2. * M_PI;

            bPhi = phi;

            return r / (2. * M_PI) * fSamp;
        }
    };
}

