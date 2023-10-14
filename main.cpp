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

std::string dec_to_string(uint64_t src) {
    std::string ans {"["};
    auto pow_exp = NEUNET_DEC_SEG_MAX / 10;
    while (pow_exp > src) {
        ans.push_back(' ');
        pow_exp /= 10;
    }
    if (src) ans += std::to_string(src);
    else ans[ans.length() - 1] = '0';
    return ans + ']';
}

uint64_t dec_div_carry(uint64_t &carry, uint64_t divd, uint64_t divr, uint64_t ans_coe) {
    auto aca = false;
    auto tmp = carry,
         ans = dec_add(aca, dec_mul(carry, ans_coe, divr), tmp);
    if (aca) ++carry;
    aca = false;
    ans = dec_sub(aca, divd, ans);
    if (aca) ++carry;
    return ans;
}
uint64_t dec_div(uint64_t divd_0, uint64_t divd_1, uint64_t divr) {
    auto     divr_d = divr / 1e19;
    uint64_t ans_0  = divd_0 / divr_d + divd_1 * 1. / divr,
             divd_2 = 0,
             divd_3 = dec_mul(divd_2, divr, ans_0);
    if (divd_0 == divd_2 && divd_1 == divd_3) return ans_0;
    auto cmp = divd_0 > divd_2 || (divd_0 == divd_2 && divd_1 > divd_3),
         tmp = false;
    if (cmp) {
        divd_3 = dec_sub(tmp, divd_1, divd_3);
        divd_2 = divd_0 - divd_2;
    } else {
        divd_3  = dec_sub(tmp, divd_3, divd_1);
        divd_2 -= divd_0;
    }
    if (tmp) --divd_2;
    uint64_t ans_1 = divd_2 / divr_d + divd_3 * 1. / divr;
    if (cmp) return ans_0 + ans_1;
    return ans_0 - ans_1 - 1;
}
uint64_t dec_div(net_set<uint64_t> &divd, const net_set<uint64_t> &divr, uint64_t &divd_high) {
    net_set<uint64_t> divd_tmp;
    if (divd.length > divr.length) {
        divd_tmp.init(divd.length - 1);
        for (auto i = 0ull; i < divr.length; ++i) {
            if (divd[i] > divr[i]) break;
            if (divd[i] < divr[i]) {
                if (divd_high && !i) break;
                divd_high = divd[0];
                divd_tmp.copy(0, divd, 1, divd_tmp.length);
                divd = std::move(divd_tmp);
                return 0;
            }
        }
    } else divd_tmp.init(divr.length);
    auto ans_coe = dec_div(divd_high, divd[0], divr[0]);
    neunet_dec_loop {
        uint64_t carry = 0;

        // if (divd_high) std::cout << dec_to_string(divd_high);
        // for (auto i = 0ull; i < divd.length; ++i) std::cout << dec_to_string(divd[i]); std::cout << '\n';
        // for (auto i = 0ull; i < divr.length; ++i) std::cout << dec_to_string(divr[i]);
        // std::cout << " * " << ans_coe << '\n';

        if (divd.length > divr.length) divd_tmp.copy(divr.length - 1, divd, divr.length, divd.length - divr.length);
        for (auto i = divr.length; i > 1; --i) {
            auto idx = i - 1;
            if (idx < divd.length) divd_tmp[idx - 1] = dec_div_carry(carry, divd[idx], divr[idx], ans_coe);
            else divd_tmp[idx - 1] = dec_div_carry(carry, 0, divr[idx], 1);
        }
        divd_high += dec_div_carry(carry, divd[0], divr[0], ans_coe);
        if (divd_high >= carry) {
            divd_high -= carry;
            break;
        }
        --ans_coe;
    }
    divd = std::move(divd_tmp);

    // if (divd_high) std::cout << dec_to_string(divd_high);
    // for (auto i = 0ull; i < divd.length; ++i) std::cout << dec_to_string(divd[i]); std::cout << '\n';
    // std::cout << std::endl;

    return ans_coe;
}
net_decimal_data dec_div_test(const net_decimal_data &divd, const net_decimal_data &divr, uint64_t prec) {
    net_assert(!dec_is_zero(divr), "neunet", "dec_div(net_decimal_data)", "Divisor could not be 0.");
    uint64_t ans_it_cnt  = 1,
             ans_idx_tmp = 0,
             divd_high   = 0;
    if (divd.it.length < divr.it.length) ans_idx_tmp = divr.it.length - divd.it.length;
    else ans_it_cnt = divd.it.length - divr.it.length + 1;
    auto ft_len = prec / NEUNET_DEC_DIG_MAX;
    if (prec % NEUNET_DEC_DIG_MAX) ++ft_len;
    net_set<uint64_t> ans_coe(ans_it_cnt + ft_len);
    auto divd_coe = divd.it,
         divr_coe = divr.it;
    divd_coe.reverse();
    divr_coe.reverse();
    divd_coe = divd_coe.unit(divd.ft);
    divr_coe = divr_coe.unit(divr.ft);
    for (auto i = ans_idx_tmp; i < ans_coe.length; ++i) ans_coe[i] = dec_div(divd_coe, divr_coe, divd_high);
    net_decimal_data ans;
    ans_idx_tmp = 0;
    while (ans_idx_tmp < ans_it_cnt && !ans_coe[ans_idx_tmp]) ++ans_idx_tmp;
    if (ans_idx_tmp < ans_it_cnt) ans.it = ans_coe.sub_set(ans_idx_tmp, ans_it_cnt - 1);
    ans.it.reverse();
    ans_idx_tmp = ans_coe.length;
    while (ans_idx_tmp > ans_it_cnt && !ans_coe[ans_idx_tmp - 1]) --ans_idx_tmp;
    if (ans_idx_tmp > ans_it_cnt) ans.ft = ans_coe.sub_set(ans_it_cnt, ans_idx_tmp - 1);
    return ans;
}

net_decimal_data dec_rem_test(net_decimal_data &divd_rem, const net_decimal_data &divr) {
    net_assert(!(divd_rem.ft.length || divr.ft.length), "neunet", "dec_rem(net_decimal_data)", "Parameters should be integer.");
    net_assert(!dec_is_zero(divr), "neunet", "dec_rem(net_decimal_data)", "Divisor could not be 0.");
    auto comp_res = dec_comp(divd_rem, divr);
    auto sgn_tmp  = false;
    if (comp_res == NEUNET_DEC_CMP_LES) return dec_init(sgn_tmp, 0);
    if (comp_res == NEUNET_DEC_CMP_EQL) return dec_init(sgn_tmp, 1);
    auto divd_ft_cnt = dec_dig_cnt(divd_rem, false),
         divd_it_cnt = dec_dig_cnt(divd_rem, true),
         divr_ft_cnt = dec_dig_cnt(divr, false),
         divr_it_cnt = dec_dig_cnt(divr, true);
    net_set<uint8_t> divr_coe(divr_it_cnt + divr_ft_cnt);
    net_set<uint8_t> divd_coe(divd_it_cnt + divd_ft_cnt);
    dec_coe(divd_coe, divd_rem);
    dec_coe(divr_coe, divr);
    net_set<uint8_t> ans_coe(divd_it_cnt - divr_it_cnt + 1);
    auto next = false;
    for (auto i = 0ull; i < ans_coe.length; ++i) ans_coe[i] = dec_div(divd_coe, divr_coe);
    divd_rem = dec_coe(divd_coe, 1);
    divd_rem.ft.reset();
    return dec_coe(ans_coe, 0);
}

int main(int argc, char *argv[], char *envp[]) {
    cout << "hello, world.\n" << endl;
    auto ch_tm_pt = NEUNET_CHRONO_TIME_POINT;

    /*
    10900471546881456|9747821456932545178
                    / 9632515584713416829
                    =   11316329001512603
    */

    auto s = false;
    auto digit = 64;
    auto a = dec_init(s, "1135147821456932545178"),
         b = dec_init(s,    "9632515584713416829"),
         c = dec_div(a, b, digit),
         d = dec_div_test(a, b, digit);
    cout << endl;
    cout << c << endl;
    cout << d << endl;
    
    uint64_t test_round = 1000000;
    auto test_time = NEUNET_CHRONO_TIME_POINT;
    for (auto i = 0ull; i < test_round; ++i) auto x = dec_div(a, b, digit);
    auto test_end = NEUNET_CHRONO_TIME_POINT;
    cout << "[normal][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;
    test_time = NEUNET_CHRONO_TIME_POINT;
    for (auto i = 0ull; i < test_round; ++i) auto x = dec_div_test(a, b, digit);
    test_end = NEUNET_CHRONO_TIME_POINT;
    cout << "[segment][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;

    /*
    117
    8454175831686503974
    4481535232212608781
    */

    // uint64_t divd  = 4316845796773021554,
    //          divr  = 9632515584713416829;
    // cout << dec_div(divd, 0, divr) << endl;

    // cout << "-0.216"_d.pow(1_d / 3_d) << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
