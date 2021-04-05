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
    struct Costas {
        static constexpr double stopBandDb = 30.0;
        Iir::Butterworth::LowPass<2> lpfI;
        Iir::Butterworth::LowPass<2> lpfQ;
        Iir::Butterworth::LowPass<2> loop;
        double dt, w, vFine, vCoarse, aI, aQ, bI, bQ;
        int lockCount{};
        int locked{false};

        Costas() = delete;

        Costas(double sampleRate, double centreFrequency, double modulationRate, double trackingRate) :
                lpfI(), lpfQ(), loop() {
            lpfI.setup(sampleRate, modulationRate);
            lpfQ.setup(sampleRate, modulationRate);
            loop.setup(sampleRate, trackingRate);
            dt = centreFrequency / sampleRate;
            w = 0.;
            vCoarse = vFine = aI = aQ = bI = bQ = 0.;
        }

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
