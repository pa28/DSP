/**
 * @file Costas.h
 * @author Richard Buckley <richard.buckley@ieee.org>
 * @version 1.0
 * @date 2021-04-04
 */

#pragma once

#include "Iir.h"

namespace rose {

    /**
     * @class Costas
     * @brief Implement a Costas Loop
     */
    struct Costas {
        Iir::Butterworth::LowPass<2> lpfI;  ///< Low pass filter for post I mixer filtering.
        Iir::Butterworth::LowPass<2> lpfQ;  ///< Low pass filter for post Q mixer filtering.
        Iir::Butterworth::LowPass<2> loop;  ///< Low pass filter for VCO control signal.

        double dt,      ///< delta time parameter for VCO free-running frequency
               w,       ///< the current phase of the VCO.
               vFine,   ///< VCO fine and fast control signal (unfiltered).
               vCoarse, ///< VCO coarse and slow control signal (filtered).
               aI,      ///< Input signal after mixing with I phase of the VCO.
               aQ,      ///< Input signal after mixing with Q phase of the VCO.
               bI,      ///< Filtered aI
               bQ;      ///< Filtered aQ

        int lockCount{};
        int locked{false};

        Costas() = delete;

        /**
         * @brief Construct and configure the Costas Loop.
         * @param sampleRate The input signal sampling rate in Hz.
         * @param centreFrequency The VCO free-running frequency in Hz.
         * @param modulationRate The modulation rate in Hz. Sets the cut-off for post mixing filters.
         * @param trackingRate The tracking rate in Hz. Sets the cut-off for the VCO control signal filter.
         */
        Costas(double sampleRate, double centreFrequency, double modulationRate, double trackingRate) :
                lpfI(), lpfQ(), loop() {
            lpfI.setup(sampleRate, modulationRate);
            lpfQ.setup(sampleRate, modulationRate);
            loop.setup(sampleRate, trackingRate);
            dt = centreFrequency / sampleRate;
            w = 0.;
            vCoarse = vFine = aI = aQ = bI = bQ = 0.;
        }

        /**
         * @brief Process an input sample.
         * @param a The input sample.
         */
        void sample(double a) {
            w = fmod(w + (2. * M_PI * (dt + vFine + vCoarse)), 2. * M_PI);
            aI = std::sin(w) * a;
            aQ = std::cos(w) * a;

            bI = lpfI.filter(aI);
            bQ = lpfQ.filter(aQ);

            vFine = bI * bQ * .06;
            auto oldCoarse = (int)(vCoarse * 100);
            vCoarse = loop.filter(vFine);
            if (oldCoarse == (int)(vCoarse * 100))
                lockCount++;
            else
                lockCount = 0;
            locked = lockCount > 10;
        }
    };
}
