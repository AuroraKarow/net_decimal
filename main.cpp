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
    
    auto sgn  = false;
    // auto divd = dec_init(sgn, "1651516516641654951811416849849811514"),
    //      divr = dec_init(sgn, "0.0000000000000000000000516511651654848");
    // dec_frac_red(divd, divr);
    // cout << divd << endl;
    // cout << divr << endl;
    
    auto num = dec_init(sgn, 241), // exp(17745) = (exp(x))^(2 * 2 * 2 * 2 * ... * 2)
         den = dec_init(sgn, 0); //,
    //      two = dec_init(sgn, 2),
    //      hlf = dec_init(sgn, 0.5),
    //      ten = dec_init(sgn, 10);
    // auto cnt = 0,
    //      bpr = 32;
    // while (dec_comp(num, ten) == NEUNET_DEC_CMP_GTR) {
    //     num = dec_mul(num, hlf);
    //     ++cnt;
    // }
    // auto ans = dec_exp(sgn, num, den, false, bpr * cnt);
    // for (auto i = 0; i < cnt; ++i) {
    //     ans = dec_mul(ans, ans);
    // }
    // dec_truncate(ans, 32);

    auto ans = dec_exp(sgn, num, den, false, 32);
    cout << ans << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
