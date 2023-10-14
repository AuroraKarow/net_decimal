/* Hello, This is Hatsune ~
 * こんにちは、ハツネちゃんです～　キラー～(∠・ω< )⌒✨
 */

#pragma once

#include <iostream>
#include "async"
#include "net_decimal"

#include <intrin.h>
#include <bitset>

#define NEUNET_CHRONO_TIME_POINT    std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()).time_since_epoch().count()

using std::cout;
using std::endl;
using std::string;

using namespace neunet;

long double __reciprocal(long double src) {
    long double curr = 0.1,
                ogn  = src;
    while (curr * ogn > 1) curr *= 0.1;
    while (curr != src) {
        src = curr;
        cout << src << endl;
        curr = src * (2 - ogn * src);
    }
    return src;
}

long double pi_iter() {
    long double base_coe = 16,
                base_dnm = 8,
                coe      = 1,
                dnm      = 0,
                ans      = 0;
    while (true) {
        ans += 1 / coe * (4 / (dnm + 1) - 2 / (dnm + 4) - 1 / (dnm + 5) - 1 / (dnm + 6));
        cout << ans << endl;
        coe *= base_coe;
        dnm += base_dnm;
    }
    return ans;
}

int main(int argc, char *argv[], char *envp[]) {
    cout << "hello, world.\n" << endl;
    auto ch_tm_pt = NEUNET_CHRONO_TIME_POINT;

    

    // 114679544494614864124 9913644494614864424
    // bool sgn = false;
    // net_decimal_base a, b;
    // a.init("114679544494614864124");
    // b.init("9913644494614864424");
    // net_decimal_frac c, d;
    // auto e = dec_add(sgn, a, true, b, true);
    // cout << e.to_string(sgn) << endl;

    // net_decimal a {-15},
    //             b {7};
    // b.modulus = true;
    // auto c = -15;
    // c %= b;
    // cout << c << endl;

    // cout << (17_d).ln() << endl;

    // net_decimal::fft_mode = true;
    // net_decimal a = -0.216,
    //             b = 1 / 3_d;
    // cout << a.dec_pow(b) << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
