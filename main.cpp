/* Hello, This is Hatsune ~
 * こんにちは、ハツネちゃんです～　キラー～(∠・ω< )⌒✨
 */

#pragma once

#define neunet_boundary_check 1

#include <iostream>
#include <fstream>
#include <bitset>
#include "net_decimal"

#define NEUNET_CHRONO_TIME_POINT    std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()).time_since_epoch().count()

using std::cout;
using std::endl;
using std::string;
using std::bitset;

using namespace neunet;

int main(int argc, char *argv[], char *envp[]) {
    cout << "hello, world.\n" << endl;
    auto ch_tm_pt = NEUNET_CHRONO_TIME_POINT;
    
    auto sgn = false;

    // auto prec = 32;

    // auto base_num = dec_init(sgn, "0.216"),
    //      base_den = dec_init(sgn, 0),
         
    //      time_num = dec_init(sgn, 1),
    //      time_den = dec_init(sgn, 3),

    //      ans = dec_series_pow(sgn, base_num, base_den, true, time_num, time_den, false, prec);

    // if (sgn) cout << '-'; cout << ans << endl;

    auto num = dec_init(sgn, "0.0001234"),
         den = dec_init(sgn, "0");

    dec_frac_red(num, den);
    cout << num << endl;
    cout << den << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
