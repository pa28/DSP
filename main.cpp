#include <iostream>
#include "fft.h"
#include "src/Costas.h"

int main() {
    const int order = 8;

    // sampling rate here 1kHz as an example
    const double sampleRate = 3000;
    const double centreFrequency = 600;
    const double modulationRate = 50;
    const double trackingRate = 10;

    rose::Costas costas{sampleRate, centreFrequency, modulationRate, trackingRate};

    double phi = 0.;
    double w = 0.;
    double dt = (centreFrequency + 20) / sampleRate;

    Iir::Butterworth::LowPass<2> lpfI;
    lpfI.setup(sampleRate, modulationRate);

    Iir::Butterworth::LowPass<2> lpfQ;
    lpfQ.setup(sampleRate, modulationRate);

    Iir::Butterworth::LowPass<2> loop;
    loop.setup(sampleRate, modulationRate);

    double lp = 0;
    for (int i = 0; i < 1000; i++) {
        if (i < 500)
            dt = (centreFrequency + 20) / sampleRate;
        else
            dt = (centreFrequency - 20) / sampleRate;
        w = fmod(w + (2.f * (double)M_PI * dt), 2. * M_PI);

        auto a = std::sin(w + phi);
        costas.sample(a);

        printf( "%4d, w%9.5f, phi%9.5f, c.w%9.5f, dw%9.5f, a%9.5f, aI%9.5f, aQ%9.5f, bI%9.5f, bQ%9.5f, vF%9.5f, vC%9.5f",
                i, w, phi, costas.w, (w+phi)-costas.w, a, costas.aI, costas.aQ, costas.bI, costas.bQ, costas.vFine, costas.vCoarse);

        std::cout << (costas.locked ? " locked\n" : "\n" );
    }

    return 0;
}
