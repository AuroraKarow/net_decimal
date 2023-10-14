/* Hello, This is Hatsune ~
 * こんにちは、ハツネちゃんです～　キラー～(∠・ω< )⌒✨
 */

#pragma once

#include <iostream>
#include "net_decimal"

#define neunet_arr_len(type, src) sizeof(src) / sizeof(type)

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

// return carry
uint64_t dec_lsh64(uint64_t &src, uint64_t bit) {
    uint64_t carry = src;
    if (bit == 0x0040) src = 0;
    else {
        carry >>= (0x0040 - bit);
        src     = uint64_t(src << bit);
    }
    auto tmp = carry;
    auto sgn = false;
    carry  <<= 1;
    src      = dec_sub(sgn, src, dec_mul(tmp, tmp, 0x158e460913d00000));
    carry   -= tmp;
    if (sgn) --carry;
    return carry;
}
void dec_lsh64(net_decimal_data &src, uint64_t bit) {
    net_decimal_data carry;
    carry.it.init(src.it.length << 1);
    auto sub_idx = 0ull;
    for (auto i = 0ull; i < src.it.length; ++i) {
        auto carry_tmp = dec_lsh64(src.it[i], bit);
        if (!carry_tmp) continue;
        sub_idx            = i + 1;
        carry.it[sub_idx] += carry_tmp % NEUNET_DEC_SEG_MAX;
        carry_tmp         /= NEUNET_DEC_SEG_MAX;
        if (carry_tmp) carry.it[++sub_idx] += carry_tmp;
    }
    if (sub_idx) carry.it = carry.it.sub_set(0, sub_idx);
    src = dec_add(carry, src);
}

void dec_lsh(net_decimal_data &src, uint64_t bit) {
    if (!(dec_is_int(src) && bit)) return;
    auto bit_cnt = bit >> 6;
    for (auto i = 0ull; i < bit_cnt; ++i) dec_lsh64(src, 0x0040);
    bit_cnt <<= 6;
    bit_cnt   = bit - bit_cnt;
    if (bit_cnt) dec_lsh64(src, bit_cnt);
}

// return remainder
net_decimal_data dec_rsh64(uint64_t &src, uint64_t bit, uint64_t src_idx, bool last) {
    auto rem_tmp = src;
    if(bit == 0x0040) src = 0;
    else {
        src   >>= bit;
        rem_tmp = (uint64_t)(rem_tmp << (0x0040 - bit));
    }
    net_decimal_data rem;
    if (last) {
        rem.ft.init(1);
        rem.ft[0] = rem_tmp;
        return rem;
    }
    if (src_idx == 1) {
        rem.it.init(1);
        rem.it[0] = rem_tmp;
        return rem;
    }
    rem.it.init(src_idx);
    rem.it[--src_idx] = rem_tmp;
    return rem;
}
void dec_rsh64(net_decimal_data &src, uint64_t bit) {
    net_decimal_data rem;
    auto src_len = src.it.length;
    for (auto i = src_len; i; --i) {
        auto idx_tmp = i - 1,
             src_tmp = src.it[idx_tmp];
        rem = dec_add(rem, dec_rsh64(src_tmp, bit, idx_tmp, i == 1));
        if (!src_tmp && src_len == i) --src_len;
        src.it[idx_tmp] = src_tmp;
    }
    auto sgn = false;
    auto b64 = dec_init(sgn, "0.542101086242752217003726400434970855712890625");
    rem      = dec_mul(rem, b64);
    if (bit == 0x0040) {
        src = std::move(rem);
        return;
    }
    if (src_len < src.it.length) src.it = src.it.sub_set(0, src_len - 1);
    src = dec_add(rem, src);
}

void dec_rsh(net_decimal_data &src, uint64_t bit) {
    if (!(dec_is_int(src) && bit)) return;
    auto bit_cnt = bit >> 6;
    for (auto i = 0ull; i < bit_cnt; ++i) dec_rsh64(src, 0x0040);
    bit_cnt <<= 6;
    bit_cnt   = bit - bit_cnt;
    if (bit_cnt) dec_rsh64(src, bit_cnt);
}

/* obsolete division */ // TODO: DIVISION

callback_dec_arg uint8_t dec_carry(int64_t &carry, const arg &coe) {
    auto tmp = carry + coe;
    carry    = 0;
    if (tmp < 0) while (tmp < 0) {
        tmp += 10;
        --carry;
    } else {
        carry = tmp / 10;
        tmp  %= 10;
    }
    return tmp;
}

void dec_coe_seg(net_set<uint64_t> &ans_coe, bool is_it, uint64_t &ans_idx, uint64_t &pow_cnt, uint64_t &seg_tmp) { if (pow_cnt == NEUNET_DEC_SEG_MAX) {
    if (is_it) ans_coe[ans_idx++] = seg_tmp;
    else ans_coe[--ans_idx] = seg_tmp;
    seg_tmp = 0;
    pow_cnt = 1;
} }

callback_dec_arg void dec_coe(net_set<arg> &dest, uint64_t &idx, const net_set<uint64_t> &seg_set, bool is_it) { for (auto i = is_it ? 0 : seg_set.length; is_it ? (i < seg_set.length) : i; is_it ? ++i : --i) {
    uint64_t seg_tmp = seg_set[is_it ? i : (i - 1)],
             seg_dig = NEUNET_DEC_DIG_MAX;
    if (i == seg_set.length) while (!(seg_tmp % 10)) {
        --seg_dig;
        seg_tmp /= 10;
    }
    while (seg_dig) {
        if (is_it && (i + 1) == seg_set.length && !seg_tmp) break;
        arg dig_tmp = 0;
        if (seg_tmp) {
            dig_tmp  = seg_tmp % 10;
            seg_tmp /= 10;
        }
        dest[--idx] = dig_tmp;
        --seg_dig;
    }
} }
// from high digit
callback_dec_arg void dec_coe(net_set<arg> &dest, const net_decimal_data &src) {
    auto idx = dest.length;
    dec_coe(dest, idx, src.ft, false);
    dec_coe(dest, idx, src.it, true);
}
// coefficient ordering from high digit
callback_dec_arg net_decimal_data dec_coe(const net_set<arg> &src, uint64_t ft_cnt) {
    auto idx_tmp = src.length;
    auto carry   = 0ll;
    while (ft_cnt && !src[idx_tmp - 1]) {
        --ft_cnt;
        --idx_tmp;
    }
    net_set<uint64_t> ans_coe(src.length);
    auto dig_cnt = ft_cnt % NEUNET_DEC_DIG_MAX,
         pow_cnt = NEUNET_DEC_SEG_MAX,
         seg_tmp = 0ull,
         ans_idx = ans_coe.length;
    if (dig_cnt) for (auto i = 0; i < dig_cnt; ++i) pow_cnt /= 10;
    else if (ft_cnt) pow_cnt = 1;
    for (auto i = 0ull; i < ft_cnt; ++i) {
        seg_tmp += dec_carry(carry, src[--idx_tmp]) * pow_cnt;
        pow_cnt *= 10;
        dec_coe_seg(ans_coe, false, ans_idx, pow_cnt, seg_tmp);
    }
    net_decimal_data ans;
    if (ft_cnt) ans.ft = ans_coe.sub_set(ans_idx, ans_coe.length - 1);
    ans_idx = 0;
    pow_cnt = 1;
    for (auto i = idx_tmp; i; --i) {
        seg_tmp += dec_carry(carry, src[i - 1]) * pow_cnt;
        pow_cnt *= 10;
        dec_coe_seg(ans_coe, true, ans_idx, pow_cnt, seg_tmp);
    }
    while (carry) {
        seg_tmp += carry % 10 * pow_cnt;
        pow_cnt *= 10;
        carry   /= 10;
        dec_coe_seg(ans_coe, true, ans_idx, pow_cnt, seg_tmp);
    }
    if (seg_tmp) ans_coe[ans_idx++] = seg_tmp;
    else while (ans_idx && !ans_coe[ans_idx - 1]) --ans_idx;
    ans.it = ans_coe.sub_set(0, ans_idx - 1);
    return ans;
}
uint8_t dec_div_obsolete(net_set<uint8_t> &divd, const net_set<uint8_t> &divr) {
    net_set<uint8_t> divd_tmp;
    if (divd.length > divr.length) {
        divd_tmp.init(divd.length - 1);
        for (auto i = 0ull; i < divr.length; ++i) {
            if (divd[i] < divr[i]) {
                divd_tmp[0] = divd[0] * 10 + divd[1];
                divd_tmp.copy(1, divd, 2, divd_tmp.length - 1);
                divd = std::move(divd_tmp);
                return 0;
            }
            if (divd[i] > divr[i]) break;
        }
    } else divd_tmp.init(divr.length);
    auto ans_coe = divd[0] / divr[0];
    if (ans_coe >= 10) 
    ans_coe = 9;
    neunet_dec_loop {
        auto carry = 0ll;
        
        // for (auto i = 0ull; i < divd.length; ++i) std::cout << dec_to_string_coe(divd[i]); std::cout << '\n';
        // for (auto i = 0ull; i < divr.length; ++i) std::cout << dec_to_string_coe(ans_coe * divr[i]); std::cout << '\n';
        
        if (divd.length > divr.length) divd_tmp.copy(divr.length - 1, divd, divr.length, divd.length - divr.length);
        for (auto i = divr.length; i > 1; --i) {
            auto idx = i - 1;
            if (idx < divd.length) divd_tmp[idx - 1] = dec_carry(carry, divd[idx] - ans_coe * divr[idx]);
            else divd_tmp[idx - 1] = dec_carry(carry, 0 - ans_coe * divr[idx]);
        }
        carry += divd[0] - ans_coe * divr[0];

        // std::cout << dec_to_string_coe(carry);
        // for (auto i = 0ull; i < divd_tmp.length; ++i) std::cout << dec_to_string_coe(divd_tmp[i]); std::cout << '\n';

        if (carry >= 0) {
            divd_tmp[0] += carry * 10;
            break;
        }

        // else std::cout << int(ans_coe) << std::endl;

        --ans_coe;
    }
    divd = std::move(divd_tmp);
    
    // std::cout << int(ans_coe) << '\n' << std::endl;
    
    return ans_coe;
}
net_decimal_data dec_div_obsolete(const net_decimal_data &divd, const net_decimal_data &divr, uint64_t prec) {
    auto divd_ft_cnt = dec_dig_cnt(divd, false),
         divd_it_cnt = dec_dig_cnt(divd, true),
         divr_ft_cnt = dec_dig_cnt(divr, false),
         divr_it_cnt = dec_dig_cnt(divr, true),
         ans_it_cnt  = 1ull,
         ans_idx_tmp = 0ull;
    net_set<uint8_t> divr_coe(divr_it_cnt + divr_ft_cnt),
                     divd_coe(divd_it_cnt + divd_ft_cnt);
    dec_coe(divd_coe, divd);
    dec_coe(divr_coe, divr);
    if (divd_it_cnt < divr_it_cnt) ans_idx_tmp = divr_it_cnt - divd_it_cnt;
    else ans_it_cnt = divd_it_cnt - divr_it_cnt + 1;
    net_set<uint8_t> ans_coe(ans_it_cnt + prec + 2);
    for (auto i = ans_idx_tmp; i < ans_coe.length; ++i) ans_coe[i] = dec_div_obsolete(divd_coe, divr_coe);

    // std::cout << "ans" << std::endl;
    // for (int tmp : ans_coe) std::cout << tmp << std::endl;

    return dec_coe(ans_coe, prec + 2);
}

net_decimal_data dec_rem_obsolete(net_decimal_data &divd_rem, const net_decimal_data &divr) {
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
    for (auto i = 0ull; i < ans_coe.length; ++i) ans_coe[i] = dec_div_obsolete(divd_coe, divr_coe);
    divd_rem = dec_coe(divd_coe, 1);
    divd_rem.ft.reset();
    return dec_coe(ans_coe, 0);
    return {};
}

// TODO: Division next of Xian Yida

uint64_t dec_div_divr(const net_decimal_data &divr, uint64_t &lsh_pow, uint64_t &seg_pow, uint64_t &sh_cnt, bool &lsh) {
    sh_cnt   = 0;
    auto ans = 0ull,
         rem = 0ull;
    auto gtr = true;
    if (divr.it.length) {        
        sh_cnt = divr.it.length - 1;
        if (sh_cnt > 1) rem = divr.it[sh_cnt - 1];
        else if (divr.ft.length) rem = divr.ft[0];
        else gtr = false;
        ans = divr.it[sh_cnt];
        lsh = false;
    } else {
        lsh = true;
        while (!divr.ft[sh_cnt++]);
        ans = divr.ft[sh_cnt - 1];
        if (sh_cnt < divr.ft.length) rem = divr.ft[sh_cnt];
        else gtr = false;
    }
    seg_pow = std::pow(10, int(std::log10(ans)) + 1);
    lsh_pow = NEUNET_DEC_SEG_MAX / seg_pow;
    if (lsh_pow > 1) {
        ans *= lsh_pow;
        ans += rem / seg_pow;
    }
    if (gtr) ++ans;
    return ans;
}

uint64_t dec_div_divd_ft_idx(const net_decimal_data &divd) {
    if (divd.it.length) return 0;
    auto ans = 0ull;
    while (!divd.ft[ans]) ++ans;
    return ans;
}

uint64_t dec_div_quot_idx_init(uint64_t divd_it_len, uint64_t sh_cnt, bool lsh, bool &quot_it) {
    auto ans = 0ull;
    if (divd_it_len) {
        quot_it = true;
        ans     = divd_it_len;
    } else {
        quot_it = false;
        ans     = 0;
    }
    if ((lsh && quot_it) || !(lsh || quot_it)) ans += sh_cnt;
    else if (ans > sh_cnt) ans -= sh_cnt;
    else {
        ans     = sh_cnt - ans;
        quot_it = !quot_it;
    }
    return ans;
}

void dec_div_quot_idx_next(uint64_t sh_cnt, uint64_t &quot_idx, bool &quot_it) {
    if (quot_it) if (quot_idx > sh_cnt) quot_idx -= sh_cnt;
    else {
        quot_it  = false;
        quot_idx = sh_cnt - quot_idx;
    } else quot_idx += sh_cnt;
}

uint64_t dec_div_divd(const net_decimal_data &divd, uint64_t &high, uint64_t lsh_pow, uint64_t seg_pow, uint64_t &it_len, uint64_t &ft_idx, uint64_t divr_seg, uint64_t &quot_idx, bool &quot_it) { neunet_dec_loop {
    auto ans = 0ull,
         rem = 0ull;
    if (it_len) {
        ans = divd.it[it_len - 1];
        if (it_len > 1) rem = divd.it[it_len - 2];
        else if (divd.ft.length) rem = divd.ft[0];
    } else if (ft_idx < divd.ft.length) {
        ans             = divd.ft[ft_idx];
        auto rem_ft_idx = ft_idx + 1;
        if (rem_ft_idx < divd.ft.length) rem = divd.ft[rem_ft_idx];
    }
    auto high_tmp = high * lsh_pow;
    high_tmp     += ans / seg_pow;
    auto ans_tmp  = ans % seg_pow;
    ans_tmp      *= lsh_pow;
    ans_tmp      += rem / seg_pow;
    if (!high_tmp && ans_tmp < divr_seg) {
        high = ans;
        dec_div_quot_idx_next(1, quot_idx, quot_it);
        if (it_len) --it_len;
        else ++ft_idx;
    }
    high = high_tmp;
    ans  = ans_tmp;
    return ans;
} }

int dec_div_next(uint64_t &ans, uint64_t &divd_0, uint64_t &divd_1, uint64_t divr) {
    ans             = divd_0 / (divr / 1e19) + divd_1 / divr * 1.;
    auto divd_0_tmp = 0ull,
         divd_1_tmp = dec_mul(divd_0_tmp, divr, ans);
    if (divd_0 == divd_0_tmp && divd_1 == divd_1_tmp) {
        divd_0 = 0;
        divd_1 = 0;
        return NEUNET_DEC_CMP_EQL;
    }
    auto sub_carry = false;
    auto divd_cmp  = 0;
    if (divd_0 > divd_0_tmp || (divd_0 == divd_0_tmp && divd_1 > divd_1_tmp)) divd_cmp = NEUNET_DEC_CMP_GTR;
    else divd_cmp = NEUNET_DEC_CMP_LES;
    if (divd_cmp == NEUNET_DEC_CMP_GTR) {
        divd_1_tmp = dec_sub(sub_carry, divd_1, divd_1_tmp);
        divd_0_tmp = divd_0 - divd_0_tmp;
    } else {
        divd_1_tmp  = dec_sub(sub_carry, divd_1_tmp, divd_1);
        divd_0_tmp -= divd_0;
    }
    if (sub_carry) --divd_0_tmp;
    divd_0 = divd_0_tmp;
    divd_1 = divd_1_tmp;
    return divd_cmp;
}
uint64_t dec_div_next(uint64_t divd_0, uint64_t divd_1, uint64_t divr /*, bool rem = false */) {
    auto ans = 0ull;
    auto cmp = NEUNET_DEC_CMP_GTR;
    neunet_dec_loop {
        auto ans_tmp = 0ull;
        if (divd_0) {
            auto cmp_tmp = dec_div_next(ans_tmp, divd_0, divd_1, divr);
            if (cmp == NEUNET_DEC_CMP_GTR) {
                ans += ans_tmp;
                cmp  = cmp_tmp;
            } else {
                ans -= ans_tmp;
                if (cmp_tmp == NEUNET_DEC_CMP_GTR) cmp = NEUNET_DEC_CMP_LES;
                else cmp = NEUNET_DEC_CMP_GTR;
            }
            if (cmp_tmp == NEUNET_DEC_CMP_EQL) return ans;
            continue;
        }
        break;
    }
    if (cmp == NEUNET_DEC_CMP_LES) {
        ans -= divd_1 / divr;
        if (divd_1 % divr) --ans;
        return ans;
    }
    return ans + divd_1 / divr;
}
net_decimal_data dec_div_next(const net_decimal_data &divd, const net_decimal_data &divr, uint64_t prec) {
    auto divd_tmp    = divd;
    auto lsh         = false,
         quot_it     = false;
    auto lsh_pow     = 0ull,
         seg_pow     = 0ull,
         sh_cnt      = 0ull,
         divr_seg    = dec_div_divr(divr, lsh_pow, seg_pow, sh_cnt, lsh),
         quot_idx    = dec_div_quot_idx_init(divd.it.length, sh_cnt, lsh, quot_it),
         divd_ft_idx = dec_div_divd_ft_idx(divd),
         divd_it_len = divd.it.length,
         divd_high   = 0ull;
    net_decimal_data ans;
    for (auto i = 0; i < prec; ++i) {
        auto divd_seg = dec_div_divd(divd_tmp, divd_high, lsh_pow, seg_pow, divd_it_len, divd_ft_idx, divr_seg, quot_idx, quot_it),
             quot_seg = dec_div_next(divd_high, divd_seg, divr_seg);
        net_decimal_data quot;
        if (quot_it) {
            quot.it.init(quot_idx--);
            quot.it[quot_idx] = quot_seg;
            if (!quot_idx) quot_it = false;
        } else {
            auto idx = quot_idx;
            quot.ft.init(++quot_idx);
            quot.ft[idx] = quot_seg;
        }
        divd_tmp = dec_add(divd_tmp, dec_mul(divr, quot), true);
        if (divd_it_len) divd_high = divd_tmp.it[--divd_it_len];
        else if (divd_ft_idx < divd_tmp.ft.length) divd_high = divd_tmp.ft[divd_ft_idx++];
        else divd_high = 0;
        if (dec_is_zero(ans)) ans = std::move(quot);
        else ans = dec_add(ans, quot);
    }
    return ans;
}

int main(int argc, char *argv[], char *envp[]) {
    cout << "hello, world.\n" << endl;
    auto ch_tm_pt = NEUNET_CHRONO_TIME_POINT;

    /*
    2222222222222222222; 16|0; 1111111111111111111|2222222222222222222
    1111111111111111112 ->
    1.9999999999999999982; 143; 9999999999999999993

    1651516516641654951811416849849811514
    0.0000000000000000000000516511651654848 ->
    31
    9744290637079135212
    0567894985125648414
    0382718827844325511.
    6887946586157766950
    2517352278691418898
    9000786681055550561
    98361501352881

    1135147821456932545178
    9632515584713416829
    ans             117.
    8454175831686503974
    mod
    8143498045462776185
    last(+)
    8143498045462776185

    11351478214569325451789632515584713416829
    69325451782515584713416829
    ans 163741857033699.
    9387570262015145428
    mod         6507975
    4955430815826696358
    last(-)      424569
    6827084768886720471
    */
    auto sgn  = false;
    auto divd = dec_init(sgn, "1135147821456932545178"),
         divr = dec_init(sgn, "9632515584713416829");
    // cout << dec_div(divd, divr, 128) << endl;
    // for (auto i = 0; i < 100000; ++i) dec_div(divd, divr, 128);
    // cout << dec_div_next(divd, divr, 11) << endl;
    // for (auto i = 0; i < 100000; ++i) dec_div_next(divd, divr, 11);

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
