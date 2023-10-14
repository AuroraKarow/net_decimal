/* Hello, This is Hatsune ~
 * こんにちは、ハツネちゃんです～　キラー～(∠・ω< )⌒✨
 */

#pragma once

#include <iostream>
#include "net_decimal"

#define NEUNET_CHRONO_TIME_POINT    std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()).time_since_epoch().count()

using std::cout;
using std::endl;
using std::string;

using namespace neunet;

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

    // net_decimal::division_precision = 64;
    // cout << ("1.35"_d).sin() << endl;

    /*
    1101001000131080019200001.00010400002170800400005100
    3179 / 256 = 12.41796875
    17 / 128 = 0.1328125
    1179 / 156 = 7.5576923076923076923076923076923
    1741376 / 128 = 13604.5
    14657816563.2 / 128 = 114514191.9
    */

    uint64_t divd_it_len = 25,
             divr_it_len = 3;
    net_set<uint8_t> divd_coe = {1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 3, 1, 0, 8, 0, 0, 1, 9, 2, 0, 0, 0, 0, 1},
                     divr_coe = {1, 5, 6};
    uint64_t dig_cnt = 1,
             seg_cnt = 300;
    if (divd_it_len >= divr_it_len) dig_cnt = divd_it_len - divr_it_len + 1;
    net_set<uint8_t> ans_set(dig_cnt + seg_cnt + 2);
    uint64_t ans_idx = 0;
    if (divd_it_len < divr_it_len) ans_idx = divr_it_len - divd_it_len;
    for (auto i = ans_idx; i < ans_set.length - 2; ++i) {
        ans_set[i] = dec_div(divd_coe, divr_coe);
        cout << '[' << (int)ans_set[i] << "][";
        for (auto tmp : divd_coe) cout << (int)tmp << ' ';
        cout << "\b]" << endl;
    }
    dig_cnt %= NEUNET_DEC_DIG_MAX;
    if (!dig_cnt) dig_cnt = NEUNET_DEC_DIG_MAX;
    uint64_t ans = 0;
    for (auto i = 0ull; i < ans_set.length - 2; ++i) {
        ans *= 10;
        ans += ans_set[i];
        if (--dig_cnt) continue;
        auto nidx = i + 1;
        while (ans_set[nidx] < 0) {
            --ans;
            ans_set[nidx] += 10;
        }
        if (!ans_set[nidx]) {
            auto carry = 0;
            while (ans_set[nidx + 1] < 0) {
                ans_set[nidx + 1] += 10;
                --carry;
            }
            while (carry < 0) {
                --ans;
                carry += 10;
            }
            ans_set[nidx] = carry;
        }
        cout << ans << '|';
        ans = 0;
        dig_cnt = NEUNET_DEC_DIG_MAX;
    }
    while (dig_cnt--) ans *= 10;
    cout << ans << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
