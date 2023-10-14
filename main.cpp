/* Hello, This is Hatsune ~
 * こんにちは、ハツネちゃんです～　キラー～(∠・ω< )⌒✨
 */

#pragma once

#include <iostream>
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
    
    // auto s = false;
    // auto a = dec_init(s, "10125.356"), // 10125|3560000000000000000
    //      b = dec_init(s,   "108.79"),  //   108|7900000000000000000
    //      c = dec_div_test(a, b, 32),
    //      d = dec_div(a, b, 32);
    // cout << c << endl;
    // cout << d << endl;
    // 93.0724882801728100009 1920213254894751355

    // net_decimal a = 17;
    // a.division_precision = 128;
    // cout << a.cos() << endl; // 102ms cos(17) = -0.2751633380515969222203426565518628
    
    // uint64_t test_round = 100;
    // auto test_time = NEUNET_CHRONO_TIME_POINT;
    // for (auto i = 0ull; i < test_round; ++i) "-0.216"_d.pow(1_d / 3_d);
    // auto test_end = NEUNET_CHRONO_TIME_POINT;
    // cout << "[test][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;
    // test_time = NEUNET_CHRONO_TIME_POINT;
    // for (auto i = 0ull; i < test_round; ++i) b = ln_4_test(32);
    // test_end = NEUNET_CHRONO_TIME_POINT;
    // cout << "[next][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;
    
    cout << "-0.216"_d.pow(1_d / 3_d) << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
