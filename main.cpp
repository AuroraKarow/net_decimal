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

/* General trigonometry series
 * num_term sine - src_num, cosine - 1
 * den_term sine - src_den, cosine - 0
 * fra_form sine - 1      , cosine - 0
 */
net_decimal_data dec_series_sin_cos(bool &ans_sgn, const net_decimal_data &src_num, const net_decimal_data &src_den, bool src_sgn, uint64_t prec, net_decimal_data &num_term, net_decimal_data &den_term, net_decimal_data &fra_form) {
    auto sers_sgn = false;
    auto one_form = dec_init(ans_sgn, 1);
    if (dec_is_zero(den_term)) den_term = one_form;
    auto num_form = dec_mul(src_num, src_num),
         den_form = dec_mul(src_den, src_den),
         fra_term = one_form,
         num      = num_term,
         den      = den_term,
         tmp      = dec_init(ans_sgn, 0);
    auto prec_seg = prec / NEUNET_DEC_DIG_MAX,
         prec_rem = prec % NEUNET_DEC_DIG_MAX;
    if (dec_is_zero(den_form)) den_form = one_form;
    ans_sgn = src_sgn;
    do {
        sers_sgn = !sers_sgn;
        num_term = dec_mul(num_term, num_form);
        den_term = dec_mul(den_term, den_form);
        for (auto i = 0; i < 2; ++i) {
            fra_form = dec_add(fra_form, one_form);
            fra_term = dec_mul(fra_term, fra_form);
        }
        tmp = dec_mul(fra_term, den_term);
        num = dec_add(ans_sgn, dec_mul(num, tmp), ans_sgn, dec_mul(num_term, den), sers_sgn != src_sgn);
        den = dec_mul(den, tmp);
    } while (dec_series_check(num_term, tmp, prec));
    return dec_div(ans_sgn, num, ans_sgn, den, false, prec);
}

int main(int argc, char *argv[], char *envp[]) {
    cout << "hello, world.\n" << endl;
    auto ch_tm_pt = NEUNET_CHRONO_TIME_POINT;
    
    auto sgn  = false;
    // auto divd = dec_init(sgn, "1651516516641654951811416849849811514"),
    //      divr = dec_init(sgn, "0.0000000000000000000000516511651654848");
    // dec_frac_red(divd, divr);
    // cout << divd << endl;
    // cout << divr << endl;
    
    auto prec = 256ull;

    auto num = dec_init(sgn, 73),
         den = dec_init(sgn, 0);

    auto fra_form = dec_init(sgn, 0),
         num_term = dec_init(sgn, 1),
         den_term = dec_init(sgn, 0);

    auto ans = dec_series_sin_cos(sgn, num, den, false, prec, num_term, den_term, fra_form);
    if (sgn) cout << '-'; cout << ans << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
