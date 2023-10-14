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

long double Aout_R(long double R) { return 205 / (1 + 20000 / R); }

struct outer {
    static struct {
        int inner_val {0};
    } inner_inst;
};

int main(int argc, char *argv[], char *envp[]) {
    cout << "hello, world.\n" << endl;
    auto ch_tm_pt = NEUNET_CHRONO_TIME_POINT;

    // 114679544494614864124.4184414 9913644494614864424.411416186156441646
    // 1135147821456932545178 9632515584713416829 117.84541758316865039744481535232
    // 114079544040094610408604012404100840414.00009900130640449040000000800644024400114106186156
    
    // auto s = false;
    // auto a = dec_init(s, "1135147821456932545178"),
    //      b = dec_init(s,    "9632515584713416829");
    // auto c = dec_rem(a, b);
    // cout << dec_to_string(false, c) << endl;
    // cout << dec_to_string(false, a) << endl;
    // auto c = dec_gcd(a, b);
    // if (c) cout << dec_to_string(false, b) << endl;
    // else cout << dec_to_string(false, a) << endl;

    // net_decimal::division_precision = 130;
    // auto pi_val = net_decimal::pi();
    // cout << pi_val << endl;

    // uint64_t test_round = 1000000;
    // auto test_time = NEUNET_CHRONO_TIME_POINT;
    // for (auto i = 0ull; i < test_round; ++i) auto c = dec_mul(a, b);
    // auto test_end = NEUNET_CHRONO_TIME_POINT;
    // cout << "[normal][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;
    // test_time = NEUNET_CHRONO_TIME_POINT;
    // for (auto i = 0ull; i < test_round; ++i) auto c = dec_fft_mul(a, b);
    // test_end = NEUNET_CHRONO_TIME_POINT;
    // cout << "[fft][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;

    // net_decimal a {-15},
    //             b {-7};
    // b.modulus = true;
    // auto c = 15;
    // c %= b;
    // cout << c << endl;

    net_decimal::division_precision = 64;
    cout << ("1.35"_d).sin() << endl;

    // net_decimal::fft_mode = true;
    // net_decimal a = -0.216,
    //             b = 1 / 3_d;
    // cout << a.dec_pow(b) << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
