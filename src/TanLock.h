/**
 * @file TanLock.h
 * @author Richard Buckley <richard.buckley@ieee.org>
 * @version 1.0
 * @date 2021-04-05
 */

#pragma once

#include "Iir.h"

namespace rose {

    /**
     * @class TanLock
     * @brief
     */
    struct TanLock {
        Iir::Butterworth::LowPass<2> lpfI;  ///< Low pass filter for post I mixer filtering.
        Iir::Butterworth::LowPass<2> lpfQ;  ///< Low pass filter for post Q mixer filtering.
        Iir::Butterworth::LowPass<2> loop;  ///< Low pass filter for VCO control signal.

        double dt,      ///< delta time parameter for VCO free-running frequency
        w,       ///< the current phase of the VCO.
        phi,     ///< The last phase output by the quadrature mixer.
        bPhi,    ///< The output of the loop control filter.
        aI,      ///< Input signal after mixing with I phase of the VCO.
        aQ,      ///< Input signal after mixing with Q phase of the VCO.
        bI,      ///< Filtered aI
        bQ;      ///< Filtered aQ

        int lockCount{};
        int locked{false};

        TanLock() = delete;

        /**
         * @brief Construct and configure the Costas Loop.
         * @param sampleRate The input signal sampling rate in Hz.
         * @param centreFrequency The VCO free-running frequency in Hz.
         * @param modulationRate The modulation rate in Hz. Sets the cut-off for post mixing filters.
         * @param trackingRate The tracking rate in Hz. Sets the cut-off for the VCO control signal filter.
         */
        TanLock(double sampleRate, double centreFrequency, double modulationRate, double trackingRate) :
                lpfI(), lpfQ(), loop() {
            lpfI.setup(sampleRate, modulationRate);
            lpfQ.setup(sampleRate, modulationRate);
            loop.setup(sampleRate, trackingRate);
            dt = centreFrequency / sampleRate;
            w = 0.;
            phi = bPhi = aI = aQ = bI = bQ = 0.;
        }

        /**
         * @brief Process an input sample.
         * @param a The input sample.
         */
        void sample(double a) {
            w = fmod(w + (2. * M_PI * (dt)) + bPhi * 0.1, 2. * M_PI);
            aI = std::sin(w) * a;
            aQ = std::cos(w) * a;

            bI = lpfI.filter(aI);
            bQ = lpfQ.filter(aQ);

            phi = std::atan2(bQ, bI);
            bPhi = loop.filter(phi);
        }
    };
}

