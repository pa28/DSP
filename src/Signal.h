/**
 * @file fft.h
 * @author Richard Buckley <richard.buckley@ieee.org>
 * @version 1.0
 * @date 2021-04-03
 */

#pragma once

#include <complex>
#include <valarray>
#include <chrono>
#include <cmath>

namespace rose {

    /**
     * @brief Constexpr function to compute the integral log base 2 of a value.
     * @param value The value to extract the log from.
     * @param l The interim log.
     * @return The log base 2 of value.
     */
    constexpr std::size_t ilog2(std::size_t value, std::size_t l = 0) {
        if (value == 1) {
            return l;
        } else if (value > 1) {
            return ilog2(value >> 1u, l + 1);
        }
        return 0;
    }

    /**
     * @brief A std::valarray based abstraction of an LTI real signal segment.
     * @tparam T The underlying data type.
     * @tparam N The number of samples in the segment.
     * @details Public inheritance from std::valarray is used for simplicity and usability, however the resize()
     * method should not be used on objects derived from rose::Signal if they change the size as this will break
     * the compile time safety of the type.
     */
    template <typename T, size_t N>
    struct Signal : public std::valarray<T> {

        std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;

        static_assert((1u << ilog2(N)) == N, "Size of Signal must be a integral power of 2.");
        static_assert(std::is_floating_point_v<T>, "Data type of Signal must be a floating point type.");

        using value_type = T;
        const std::size_t size{N};

        constexpr Signal() noexcept : std::valarray<T>(N) {}

        /**
         * @brief Fill the signal segment with a sinusoidal wave form of a given frequency, amplitude and phase.
         * @details The frequency is the number of cycles completed in the length of the segment.
         * @param f The signal frequency comuted by the desired frequency divided by the sample rate.
         * @param amp The signal amplitude.
         * @param phi The signal phase in radians.
         */
        void sinusoid(T f, T amp = 1., T phi = 0.) {
            T t = 0;
            for (auto &s : *this) {
                T w = (2. * M_PI * t) + phi;
                s = amp * sin(w);
                t += f;
            }
        }

        /**
         * @brief Apply a Hann window to the signal.
         */
        void hann() {
            T h = 0.;
            T n = 0.;
            for (auto &s : *this) {
                h = sin((M_PI * n)/size);
                s *= h * h;
                n += 1;
            }
        }
    };

    /**
     * @brief A std::valarray based abstraction of an LTI complex signal segment.
     * @tparam T The underlying data type.
     * @tparam N The number of samples in the segment.
     * @details Public inheritance from std::valarray is used for simplicity and usability, however the resize()
     * method should not be used on objects derived from rose::ComplexSignal if they change the size as this will break
     * the compile time safety of the type.
     */
    template <typename T, size_t N>
    struct ComplexSignal : public std::valarray<std::complex<T>> {
    protected:
        std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;

        // Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
        // Better optimized but less intuitive
        // Derived from https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B
        // !!! Warning : in some cases this code make result different from not optimased version above (need to fix bug)
        // The bug is now fixed @2017/05/30

        /**
         * @brief Compute the complex Fast Fourier Transform of a std::valarray<std::complex<T>>.
         * @details The transform is computed in place and is not normalized.
         * @param x The std::valarray.
         */
        void fft_dif(std::valarray<std::complex<T>>& x)
        {
            // DFT
            unsigned int k = N, n;
            double thetaT = 3.14159265358979323846264338328L / N;
            Complex phiT = Complex(std::cos(thetaT), -std::sin(thetaT)), T0;
            while (k > 1)
            {
                n = k;
                k >>= 1u;
                phiT = phiT * phiT;
                T0 = 1.0L;
                for (unsigned int l = 0; l < k; l++)
                {
                    for (unsigned int a = l; a < N; a += n)
                    {
                        unsigned int b = a + k;
                        Complex t = x[a] - x[b];
                        x[a] += x[b];
                        x[b] = t * T0;
                    }
                    T0 *= phiT;
                }
            }
            // Decimate
            auto m = (unsigned int)log2(N);
            for (unsigned int a = 0; a < N; a++)
            {
                unsigned int b = a;
                // Reverse bits
                b = (((b & 0xaaaaaaaau) >> 1u) | ((b & 0x55555555u) << 1u));
                b = (((b & 0xccccccccu) >> 2u) | ((b & 0x33333333u) << 2u));
                b = (((b & 0xf0f0f0f0u) >> 4u) | ((b & 0x0f0f0f0fu) << 4u));
                b = (((b & 0xff00ff00u) >> 8u) | ((b & 0x00ff00ffu) << 8u));
                b = ((b >> 16u) | (b << 16u)) >> (32u - m);
                if (b > a)
                {
                    Complex t = x[a];
                    x[a] = x[b];
                    x[b] = t;
                }
            }
            //// Normalize (This section make it not working correctly)
            //Complex f = 1.0 / sqrt(N);
            //for (unsigned int i = 0; i < N; i++)
            //	x[i] *= f;
            stop = std::chrono::high_resolution_clock::now();
        }

        /**
         * @brief Compute the inverse complex Fast Fourier Transform.
         * @details See fft_dft()
         */
        void ifft_dif(std::valarray<std::complex<T>>& x)
        {
            // conjugate the complex numbers
            x = x.apply(std::conj);

            // forward fft
            fft_dif( x );

            // conjugate the complex numbers again
            x = x.apply(std::conj);

            // scale the numbers
            x /= x.size();
        }


    public:
        static_assert((1u << ilog2(N)) == N, "Size of Signal must be a integral power of 2.");
        static_assert(std::is_floating_point_v<T>, "Data type of Signal must be a floating point type.");

        using value_type = T;
        const std::size_t size{N};
        using Complex = std::complex<T>;

        /**
         * @brief Construct a complex signal segment with N complex samples.
         * @details N must be an integral power of 2
         */
        constexpr ComplexSignal() noexcept : std::valarray<std::complex<T>>(N) {}

        /**
         * @brief Assign a real signal segment to a this complex segment.
         * @details The real segment must have N*2 samples. Even samples are assigned to the real values, add
         * samples are assigned to imaginary values. This assignment allows the FFT N*2 long real signal segment
         * to be computed by a complex FFT of N complex samples.
         * @param real The real segment.
         * @return A reference to this complex segment.
         */
        constexpr ComplexSignal<T,N>& operator=(const Signal<T,N*2> &real) noexcept {
            for (std::size_t i = 0; i < N; ++i) {
                (*this)[i] = std::complex<T>{real[i*2], real[i*2+1]};
            }
            return *this;
        }

        /**
         * @brief Fill the signal segment with a complex sinusoidal wave form of a given frequency, amplitude and phase.
         * @details The frequency is the number of cycles completed in the length of the segment.<p/>
         * The start and stop times of the function are captured using the system high resolution clock. The
         * duration of the function call can be found by calling duration() before the next call to an instrumentd
         * function.
         * @param f The signal frequency.
         * @param amp The signal amplitude.
         * @param phi The signal phase in radians.
         */
        void sinusoid(T f, T amp = 1.0, T phi = 0.) {
            T t = 0;
            start = std::chrono::high_resolution_clock::now();
            for (auto &s : *this) {
                T w = (2. * M_PI * f * t) + phi;
                s = Complex{amp * std::cos(w), amp * std::sin(w)};
                t += 1/(T)N;
            }
            stop = std::chrono::high_resolution_clock::now();
        }

        /**
         * @brief Apply a Hann window to the signal.
         */
        void hann() {
            T h = 0.;
            T n = 0.;
            for (auto &s : *this) {
                h = sin((M_PI * n)/size);
                s *= h * h;
                n += 1;
            }
        }

        /**
         * @brief Compute an in-place complex Fast Fourier Transform of the signal.
         * @details No normalization is perform.<p/>
         * The start and stop times of the function are captured using the system high resolution clock. The
         * duration of the function call can be found by calling duration() before the next call to an instrumentd
         * function.
         */
        void fft() {
            start = std::chrono::high_resolution_clock::now();
            fft_dif(*this);
            stop = std::chrono::high_resolution_clock::now();
        }

        /**
         * @brief Compute an in-place inverse complex Fast Fourier Transform of the signal.
         * @details
         * The start and stop times of the function are captured using the system high resolution clock. The
         * duration of the function call can be found by calling duration() before the next call to an instrumentd
         * function.
         */
        void ifft() {
            start = std::chrono::high_resolution_clock::now();
            ifft_dif(*this);
            stop = std::chrono::high_resolution_clock::now();
        }

        /**
         * @brief Compute an return the duration of the last instrumented function.
         * @return std::chrono::duration_cast<std::chrono::microseconds>
         */
        auto duration() const {
            return std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        }
    };
}

