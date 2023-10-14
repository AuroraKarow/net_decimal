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

// net_decimal ln_proto(uint64_t prec, bool show_iter = false) {
//     net_decimal b {1},
//                 p {0},
//                 q {b},
//                 o {p},
//                 u {*this - 1},
//                 v {*this + 1},
//                 h {u * u},
//                 s {v * v},
//                 m {p},
//                 n {p};
//     do {
//         if (show_iter) std::cout << n << std::endl;
//         m = std::move(n);
//         if (q == v) p += u;
//         else {
//             o  = b * v;
//             p  = o * p + u * q;
//             q *= o;
//         }
//         u *= h;
//         v *= s;
//         b += 2;
//         n.base = dec_div(n.sgn, p.base, p.sgn, q.base, q.sgn, prec);
//     } while (m != n);
//     return 2 * n;
// }

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
    return ans_0 - ans_1;
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

    // std::cout << '\n';
    // if (divd_high) std::cout << dec_to_string_seg(divd_high);
    // for (auto i = 0ull; i < divd.length; ++i) std::cout << dec_to_string_seg(divd[i]); std::cout << '\n';
    // std::cout << std::endl;
    
    for (auto i = 0ull; ; ++i) {
        uint64_t carry = 0;
        if (divd.length > divr.length) divd_tmp.copy(divr.length - 1, divd, divr.length, divd.length - divr.length);
        for (auto i = divr.length; i > 1; --i) {
            auto idx = i - 1;
            if (idx < divd.length) divd_tmp[idx - 1] = dec_div_carry(carry, divd[idx], divr[idx], ans_coe);
            else divd_tmp[idx - 1] = dec_div_carry(carry, 0, divr[idx], 1);
        }
        auto tmp = dec_div_carry(carry, divd[0], divr[0], ans_coe);

        // if (tmp) std::cout << dec_to_string_seg(tmp);
        // for (auto i = 0ull; i < divd_tmp.length; ++i) std::cout << dec_to_string_seg(divd_tmp[i]);
        // std::cout << ans_coe;
        
        if (carry <= divd_high) {
            if (ans_pow == 1 || !i) {
                divd_high = tmp;
                break;
            }
            ans_coe += ans_pow;
            ans_pow /= 10;
        }
        
        // else std::cout << '*';
        // std::cout << '\n';

        ans_coe -= ans_pow;
    }
    divd = std::move(divd_tmp);

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

net_decimal_data dec_e10_mul(const net_decimal_data &src, uint64_t dig_cnt) {
    auto seg_cnt = dig_cnt / NEUNET_DEC_DIG_MAX,
         seg_rem = 0ull;
    net_decimal_data ans;
    if (seg_cnt) {
        if (seg_cnt >= src.ft.length) if (seg_cnt == src.ft.length) ans.it = std::move(src.ft);
        else {
            ans.it.init(seg_cnt);
            ans.it.copy(0, src.ft, 0, src.ft.length);
        } else {
            ans.ft = src.ft.sub_set(seg_cnt, src.ft.length - 1);
            ans.it = src.ft.sub_set(0, seg_cnt - 1);
        }
        ans.it.reverse();
        ans.it = ans.it.unit(src.it);
    }
    dig_cnt %= NEUNET_DEC_DIG_MAX;
    if (!dig_cnt) return ans;
    seg_cnt = 1;
    for (auto i = 0; i < dig_cnt; ++i) seg_cnt *= 10;
    dig_cnt = NEUNET_DEC_SEG_MAX / seg_cnt;
    for (auto i = ans.ft.length; i; --i) {
        auto seg_idx     = i - 1;
        auto seg_tmp     = ans.ft[seg_idx] / dig_cnt;
        ans.ft[seg_idx] %= dig_cnt;
        ans.ft[seg_idx] *= seg_cnt;
        ans.ft[seg_idx] += seg_rem;
        seg_rem          = seg_tmp;
    }
    if (ans.ft.length) {
        auto tmp = ans.ft.length;
        while (tmp && !ans.ft[tmp - 1]) --tmp;
        if (tmp) ans.ft = ans.ft.sub_set(0, tmp - 1);
    }
    for (auto i = 0ull; i < ans.it.length; ++i) {
        auto seg_tmp = ans.it[i] / dig_cnt;
        ans.it[i]   %= dig_cnt;
        ans.it[i]   *= seg_cnt;
        ans.it[i]   += seg_rem;
        seg_rem      = seg_tmp;
    }
    if (seg_rem) {
        net_set<uint64_t> it(ans.it.length + 1);
        it.copy(0, ans.it, 0, ans.it.length);
        it[ans.it.length] = seg_rem;
        ans.it = std::move(it);
    }
    return ans;
}

net_decimal_data dec_e10_div(const net_decimal_data &src, uint64_t dig_cnt, bool div = false) {
    auto seg_cnt = dig_cnt / NEUNET_DEC_DIG_MAX,
         seg_rem = 0ull;
    net_decimal_data ans;       
    if (seg_cnt) {
        if (seg_cnt >= src.it.length) if (seg_cnt == src.it.length) ans.ft = std::move(src.it);
        else {
            ans.ft.init(seg_cnt);
            ans.ft.copy(0, src.it, 0, src.it.length);
        } else {
            ans.it = src.it.sub_set(seg_cnt, src.it.length - 1);
            ans.ft = src.it.sub_set(0, seg_cnt - 1);
        }
        ans.ft.reverse();
        ans.ft = ans.ft.unit(src.ft);
    }
    dig_cnt %= NEUNET_DEC_DIG_MAX;
    if (!dig_cnt) return ans;
    seg_cnt = 1;
    for (auto i = 0; i < dig_cnt; ++i) seg_cnt *= 10;
    dig_cnt = NEUNET_DEC_SEG_MAX / seg_cnt;
    for (auto i = ans.it.length; i; --i) {
        auto seg_idx     = i - 1;
        auto seg_tmp     = ans.it[seg_idx] % seg_cnt;
        ans.it[seg_idx] /= seg_cnt;
        ans.it[seg_idx] += seg_rem * dig_cnt;
        seg_rem          = seg_tmp;
    }
    if (ans.it.length) {
        auto tmp = ans.it.length;
        while (tmp && !ans.it[tmp - 1]) --tmp;
        if (tmp) ans.it = ans.it.sub_set(0, tmp - 1);
    }
    for (auto i = 0ull; i < ans.ft.length; ++i) {
        auto seg_tmp = ans.ft[i] % seg_cnt;
        ans.ft[i]   /= seg_cnt;
        ans.ft[i]   += seg_rem * dig_cnt;
        seg_rem      = seg_tmp;
    }
    if (seg_rem) {
        net_set<uint64_t> ft(ans.ft.length + 1);
        ft.copy(0, ans.ft, 0, ans.ft.length);
        ft[ans.ft.length] = seg_rem;
        ans.ft = std::move(ft);
    }
    return ans;
}

int main(int argc, char *argv[], char *envp[]) {
    cout << "hello, world.\n" << endl;
    auto ch_tm_pt = NEUNET_CHRONO_TIME_POINT;

    // 114679544494614864124.4184414 9913644494614864424.411416186156441646
    // 1135147821456932545178 9632515584713416829 117.84541758316865039744481535232
    // 114079544040094610408604012404100840414.00009900130640449040000000800644024400114106186156
    
    auto s = false;
    auto a = dec_init(s, "1135147821456932545178"),
         b = dec_init(s,    "9632515584713416829");
    auto c = dec_rem(a, b);
    cout << c << endl;
    cout << a << endl;

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

    // uint64_t divd_it_len = 25,
    //          divr_it_len = 3;
    // net_set<uint8_t> divd_coe = {1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 3, 1, 0, 8, 0, 0, 1, 9, 2, 0, 0, 0, 0, 1},
    //                  divr_coe = {1, 5, 6};
    // uint64_t dig_cnt = 1,
    //          seg_cnt = 300;
    // if (divd_it_len >= divr_it_len) dig_cnt = divd_it_len - divr_it_len + 1;
    // net_set<uint8_t> ans_set(dig_cnt + seg_cnt + 2);
    // uint64_t ans_idx = 0;
    // if (divd_it_len < divr_it_len) ans_idx = divr_it_len - divd_it_len;
    // for (auto i = ans_idx; i < ans_set.length - 2; ++i) {
    //     ans_set[i] = dec_div(divd_coe, divr_coe);
    //     cout << '[' << (int)ans_set[i] << "][";
    //     for (auto tmp : divd_coe) cout << (int)tmp << ' ';
    //     cout << "\b]" << endl;
    // }
    // dig_cnt %= NEUNET_DEC_DIG_MAX;
    // if (!dig_cnt) dig_cnt = NEUNET_DEC_DIG_MAX;
    // uint64_t ans = 0;
    // for (auto i = 0ull; i < ans_set.length - 2; ++i) {
    //     ans *= 10;
    //     ans += ans_set[i];
    //     if (--dig_cnt) continue;
    //     auto nidx = i + 1;
    //     while (ans_set[nidx] < 0) {
    //         --ans;
    //         ans_set[nidx] += 10;
    //     }
    //     if (!ans_set[nidx]) {
    //         auto carry = 0;
    //         while (ans_set[nidx + 1] < 0) {
    //             ans_set[nidx + 1] += 10;
    //             --carry;
    //         }
    //         while (carry < 0) {
    //             --ans;
    //             carry += 10;
    //         }
    //         ans_set[nidx] = carry;
    //     }
    //     cout << ans << '|';
    //     ans = 0;
    //     dig_cnt = NEUNET_DEC_DIG_MAX;
    // }
    // while (dig_cnt--) ans *= 10;
    // cout << ans << endl;

    // auto sgn  = false;
    // auto test = dec_init(sgn, "324938271560444852566064687276808488928"),
    //      four = dec_init(sgn, 4);
    // cout << dec_div_test(test, four, 0) << endl;

    // uint64_t test_round = 10000000;
    // auto test_time = NEUNET_CHRONO_TIME_POINT;
    // for (auto i = 0ull; i < test_round; ++i) dec_div_test(test, four, 0);
    // auto test_end = NEUNET_CHRONO_TIME_POINT;
    // cout << "[normal][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
