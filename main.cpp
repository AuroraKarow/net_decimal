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

// true - divdend >= next dividend; false - dividend < next divdend
bool dec_div_ge(uint64_t divd_high, const net_set<uint64_t> &divd, uint64_t divd_tmp_carry, uint64_t divd_high_tmp, const net_set<uint64_t> &divd_tmp) {
    if (divd_high > divd_tmp_carry) return true;
    if (divd_high < divd_tmp_carry) return false;
    if (divd[0] > divd_high_tmp) return true;
    if (divd[0] < divd_high_tmp) return false;
    for (auto i = 1ull; i < divd.length; ++i) {
        auto divd_tmp_coe = divd_tmp[i - 1];
        if (divd[i] > divd_tmp_coe) return true;
        if (divd[i] < divd_tmp_coe) return false;
    }
    if (divd_tmp[divd_tmp.length - 1]) return false;
    return true;
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
    /*
    4316845796773021554|0000000000000000000
                      / 9632515584713416829
                      = 4481535232212608781
    precision lost -> 4316845796773021554 / 0.9632515584713416829 + 0
                      = 4481535232212609024
                      * 9632515584713416829
    4316845796773021787|9395928199718864896
    or return 4481535232212609024
    */
    auto     divr_d = divr / 1e19;
    uint64_t ans_0  = divd_0 / divr_d + divd_1 * 1. / divr,
             divd_2 = 0,
             divd_3 = dec_mul(divd_2, divr, ans_0);
    if (divd_0 == divd_2 && divd_1 == divd_3) return ans_0;
    /*
    compare = 4316845796773021554|0000000000000000000 > 4316845796773021787|9395928199718864896
    true  4316845796773021554|0000000000000000000 - 4316845796773021787|9395928199718864896
    false 4316845796773021787|9395928199718864896 - 4316845796773021554|0000000000000000000
    = 233|9395928199718864896
    */
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
    /*
    242 = 233 / 0.9632515584713416829 + 9395928199718864896 / 9632515584713416829
    compare ->
    true  4481535232212609024 + 242
    truncating for rounding down
    -> false 4481535232212609024 - (242 + 1)
    Rounding up
    return 4481535232212608781
    */
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
    auto ans_coe = dec_div(divd_high, divd[0], divr[0]),
         ans_pow = NEUNET_DEC_SEG_MAX;
    while (ans_coe < ans_pow) ans_pow /= 10;

    std::cout << '\n';
    if (divd_high) std::cout << dec_to_string_seg(divd_high);
    for (auto i = 0ull; i < divd.length; ++i) std::cout << dec_to_string_seg(divd[i]); std::cout << '\n';
    std::cout << std::endl;
        
    for (auto i = 0ull; ; ++i) {
        uint64_t carry = 0;
        if (divd.length > divr.length) divd_tmp.copy(divr.length - 1, divd, divr.length, divd.length - divr.length);
        for (auto i = divr.length; i > 1; --i) {
            auto idx = i - 1;
            if (idx < divd.length) divd_tmp[idx - 1] = dec_div_carry(carry, divd[idx], divr[idx], ans_coe);
            else divd_tmp[idx - 1] = dec_div_carry(carry, 0, divr[idx], 1);
        }
        auto tmp = dec_div_carry(carry, divd[0], divr[0], ans_coe);

        if (tmp) std::cout << dec_to_string_seg(tmp);
        for (auto i = 0ull; i < divd_tmp.length; ++i) std::cout << dec_to_string_seg(divd_tmp[i]);
        std::cout << ans_coe;
        
        if (carry <= divd_high) {
            if (ans_pow == 1 || !i) {
                divd_high = tmp;
                break;
            }
            ans_coe += ans_pow;
            ans_pow /= 10;
        }
        
        else std::cout << '*';
        std::cout << '\n';

        ans_coe -= ans_pow;
    }
    divd = std::move(divd_tmp);

    std::cout << std::endl;

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

int main(int argc, char *argv[], char *envp[]) {
    cout << "hello, world.\n" << endl;
    auto ch_tm_pt = NEUNET_CHRONO_TIME_POINT;

    /*
    10900471546881456|9747821456932545178
                    / 9632515584713416829
                    =   11316329001512603
    */

    // auto s = false;
    // auto digit = 64;
    // auto a = dec_init(s, "1135147821456932545178"),
    //      b = dec_init(s,    "9632515584713416829"),
    //      c = dec_div(a, b, digit),
    //      d = dec_div_test(a, b, digit);
    // cout << endl;
    // cout << c << endl;
    // cout << d << endl;

    /*
    109004715468814569747821456932545178748411534634
                         / 1646417416466461649632515
    a = 109004715468814569747821456932545178748411534634 / 164641
    b = 109004715468814569747821456932545178748411534634 / 7416466461649632515
    return a * NEUNET_DEC_SEG_MAX + b
    */

    auto s = false;
    auto a = dec_init(s, "10125.356"),
         b = dec_init(s,   "108.79"),
         c = dec_div_test(a, b, 32),
         d = dec_div(a, b, 32);
    cout << c << endl;
    cout << d << endl;

    // uint64_t test_round = 1000000;
    // auto test_time = NEUNET_CHRONO_TIME_POINT;
    // for (auto i = 0ull; i < test_round; ++i) c = dec_div(a, b, 32);
    // auto test_end = NEUNET_CHRONO_TIME_POINT;
    // cout << "[normal][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;
    // test_time = NEUNET_CHRONO_TIME_POINT;
    // for (auto i = 0ull; i < test_round; ++i) auto c = dec_div_test(a, b, 32);
    // test_end = NEUNET_CHRONO_TIME_POINT;
    // cout << "[segment][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;

    /*
    7 8860000000000000000|0000000000000000000 / 108|7900000000000000000
                         = 724882801728100009
    ----------------------
       730185185185185185 --ans -> 724882801728100009
       630185185185185185
    */ 
    

    // cout << "-0.216"_d.pow(1_d / 3_d) << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
