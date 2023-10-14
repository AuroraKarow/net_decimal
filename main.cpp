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

    // auto sgn  = false;
    // auto test = dec_init(sgn, "84344328648949415486748944871499497167"),
    //      four = dec_init(sgn, 4),
    //      ans  = dec_mul(test, four);
    // cout << ans << endl;
    
    // cout << "-0.216"_d.pow(1_d / 3_d) << endl;

    // uint64_t test_round = 10000;
    // auto test_time = NEUNET_CHRONO_TIME_POINT;
    // for (auto i = 0ull; i < test_round; ++i) "-0.216"_d.pow(1_d / 3_d);
    // auto test_end = NEUNET_CHRONO_TIME_POINT;
    // cout << "[test][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;

    cout << net_decimal_ln4.get(128) << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
