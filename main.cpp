/* Hello, This is Hatsune ~
 * こんにちは、ハツネちゃんです～　キラー～(∠・ω< )⌒✨
 */

#pragma once

#define neunet_boundary_check 1

#include <iostream>
#include <fstream>
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
    // auto divd = dec_init(sgn, "1651516516641654951811416849849811514"),
    //      divr = dec_init(sgn, "0.0000000000000000000000516511651654848");
    // dec_frac_red(divd, divr);
    // cout << divd << endl;
    // cout << divr << endl;

    auto ans = 176.145_d / 113.2_d;
    cout << double(ans) << endl;

    auto fst = 1445_d / 85_d, snd = 22_d;
    cout << (fst | snd) << endl;
    cout << (17 | 22) << endl;
    cout << fst.absolute << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
