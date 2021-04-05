#include <iostream>
#include "Costas.h"
#include "fft.h"

int main() {
    const int order = 8;

    // sampling rate here 1kHz as an example
    const float sampleRate = 3000;
    const float centreFrequency = 200;
    const float modulationRate = 50;
    const float bandStopDb = 30;

    rose::Costas<double> costas{sampleRate, centreFrequency, modulationRate, bandStopDb, modulationRate, bandStopDb};

    float phi = M_PI_4;
    float w = 0.;
    float wC = 0.;
    float dt = (centreFrequency - 1) / sampleRate;
    float dtC = (centreFrequency) / sampleRate;

    Iir::Butterworth::LowPass<2> lpfI;
    lpfI.setup(sampleRate, modulationRate);

    Iir::Butterworth::LowPass<2> lpfQ;
    lpfQ.setup(sampleRate, modulationRate);

    Iir::Butterworth::LowPass<2> loop;
    loop.setup(sampleRate, modulationRate);

    float lp = 0;
    for (int i = 0; i < 1000; i++) {
        w = w + (2.f * (float)M_PI * dt);
        wC = wC + (2.f * (float)M_PI * (dtC + lp));

        auto a = std::sin(w );

        auto aI = std::sin(wC) * a;
        auto aQ = std::cos(wC) * a;

        auto bI = lpfI.filter(aI);
        auto bQ = lpfQ.filter(aQ);

        lp = loop.filter(bI * bQ) / 5;
        lp = bI * bQ / 18;

        printf( "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n", w, wC, a, aI, aQ, bI, bQ, bI * bQ / 18, lp);

        w = fmod(w, 2. * M_PI);
        wC = fmod(wC, 2. * M_PI);
    }

    return 0;
}
