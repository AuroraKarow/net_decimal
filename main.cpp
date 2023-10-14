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

    // 1135147821456932545178.9632515584713416829 117.84541758316865039744481535232

    // auto sgn  = false;
    // auto test = dec_init(sgn, "84344328648949415486748944871499497167"),
    //      four = dec_init(sgn, 4),
    //      ans  = dec_mul(test, four);
    // cout << ans << endl;
    
    // cout << std::pow(-0.216_d, (1_d / 3_d)) << endl;

    // uint64_t test_round = 10000;
    // auto test_time = NEUNET_CHRONO_TIME_POINT;
    // for (auto i = 0ull; i < test_round; ++i) (17_d).ln();
    // auto test_end = NEUNET_CHRONO_TIME_POINT;
    // cout << "[test][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;

    auto test = 42163549464864328_d;
    cout << test.cos() << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
