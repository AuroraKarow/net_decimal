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

    net_decimal::default_division_precision = 32;

    cout << ("192831284698325983253481263856325632"_d >> 64_d) << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
