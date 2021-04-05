
#define COSTAS 0
#define TANLOCK 0
#define FDESC 1


#include <iostream>
#include <iomanip>
#if COSTAS
#include "src/Costas.h"
#endif
#if TANLOCK
#include "src/TanLock.h"
#endif
#if FDESC
#include "src/FreqDesc.h"
#endif

int main() {
    const int order = 8;

    // sampling rate here 1kHz as an example
    const double sampleRate = 3000;
    const double centreFrequency = 600;
    const double modulationRate = 50;
    const double trackingRate = 10;

#if COSTAS
    rose::Costas costas{sampleRate, centreFrequency, modulationRate, trackingRate};

    double phi = 0.;
    double w = 0.;
    double dt = (centreFrequency + 20) / sampleRate;

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

#endif
#if TANLOCK

    rose::TanLock tanLock{sampleRate, centreFrequency, modulationRate, modulationRate * 2};
    double phi = 0.;
    double w = 0.;
    double dt = (centreFrequency + 20) / sampleRate;

    for (int i = 0; i < 1000; i++) {
        if (i < 500)
            dt = (centreFrequency + 50) / sampleRate;
        else
            dt = (centreFrequency - 50) / sampleRate;
        w = fmod(w + (2.f * (double) M_PI * dt), 2. * M_PI);

        auto a = std::sin(w + phi);
        tanLock.sample(a);

        printf( "%4d, w%9.5f, phi%9.5f, c.w%9.5f, dw%9.5f, a%9.5f, aI%9.5f, aQ%9.5f, bI%9.5f, bQ%9.5f, tl.phi%9.5f, tl.bPhi%9.5f\n",
                i, w, phi, tanLock.w, (w+phi)-tanLock.w, a, tanLock.aI, tanLock.aQ, tanLock.bI, tanLock.bQ, tanLock.phi, tanLock.bPhi);
//        std::cout << (2. * M_PI * (dt)) << ' ' << (2. * M_PI * (tanLock.dt)) + tanLock.phi * 0.1 << '\n';
    }

#endif
#if FDESC
    rose::FreqDesc freqDesc{sampleRate, centreFrequency, modulationRate};
    double phi = 0.;
    double w = 0.;
    double dt = 0.;

    for (int i = 0; i < 1000; i++) {
        if (i < 500)
            dt = (centreFrequency + 50) / sampleRate;
        else
            dt = (centreFrequency - 50) / sampleRate;
        w = std::fmod(w + (2. * M_PI * (dt)), 2. * M_PI);

        auto a = std::sin(w + phi);
        auto b = freqDesc.sample(a);

        printf( "%4d, w%9.5f, phi%9.5f, c.w%9.5f, a%9.5f, b%9.5f, aI%9.5f, aQ%9.5f, bI%9.5f, bQ%9.5f, tl.phi%9.5f, tl.bPhi%9.5f\n",
                i, w, phi, freqDesc.w, a, b, freqDesc.aI, freqDesc.aQ, freqDesc.bI, freqDesc.bQ, freqDesc.phi, freqDesc.bPhi);
    }

#endif

    return 0;
}
