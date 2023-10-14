/* Hello, This is Hatsune ~
 * こんにちは、ハツネちゃんです～　キラー～(∠・ω< )⌒✨
 */

#pragma once

#include <iostream>
#include <fstream>
#include "net_decimal"

#define NEUNET_CHRONO_TIME_POINT    std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()).time_since_epoch().count()

using std::cout;
using std::endl;
using std::string;

using namespace neunet;

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

#define neunet_arr_len(type, src) sizeof(src) / sizeof(type)

int main(int argc, char *argv[], char *envp[]) {
    cout << "hello, world.\n" << endl;
    auto ch_tm_pt = NEUNET_CHRONO_TIME_POINT;

    // 106919138837687771558617554786784057736000018581285323637421131169792
    // >> 100 0.0000000000000000000000000000007888609052210118054117285652827862296732064351090230047702789306640625
    // << 100 1267650600228229401496703205376
    // >> 64  0.0000000000000000000542101086242752217003726400434970855712890625
    // << 64  18446744073709551616

    auto sgn  = false;
    auto test = dec_init(sgn, "84344328648949415486748944871499497167"),
         lshb = dec_init(sgn, "1267650600228229401496703205376"),
         rshb = dec_init(sgn, 
         "0.0000000000000000000000000000007888609052210118054117285652827862296732064351090230047702789306640625");

    /*
    106919138837 6877715586175547867 8405773600001858128 5323637421131169792 >> 64
                         57960981304 0501766838244223336 9556205003354406912

                         57960981303 6773349748109751544 6890611201524734497.0703125
                                     3728417090134471791 8108814802558015344.4025199860334396362 3046875
                                                         4556778999271657070.2385725510907832358 4981262683868408203 125
                                                                           0.2885949628757771279 1971862316131591796 875
    */

    cout << test << endl;

    test = dec_mul(test, lshb);
    test = dec_mul(test, rshb);
    cout << test << endl;

    dec_lsh(test, 100);
    dec_rsh(test, 100);
    cout << test << endl;

    // uint64_t test_round = 1000000;
    // auto test_time = NEUNET_CHRONO_TIME_POINT;
    // for (auto i = 0ull; i < test_round; ++i) {
    //     dec_lsh(test, 100);
    //     dec_rsh(test, 100);
    // }
    // auto test_end = NEUNET_CHRONO_TIME_POINT;
    // cout << "[mul][" << test_round << "][" << test_end - test_time << "ms]" << std::endl;

    // 1135147821456932545178.9632515584713416829 117.84541758316865039744481535232
    
    // cout << std::pow(-0.216_d, (1_d / 3_d)) << endl;

    // auto test = "1022183956076299059836405850183740184811657306598078"_d / "12119178283234"_d;
    // auto test = "0.216"_d;
    // // test.division_precision = 128;
    // cout << test.exp() << endl;

    // auto src = "84344328648949415486748944871499497167"_d,
    //      fra = 0.5_d;
    // auto cnt = 0ull;
    // while (src > 10) {
    //     src *= fra;
    //     ++cnt;
    // }
    // auto src_prec = src.division_precision,
    //      src_fold = (uint64_t)std::log2(cnt);
    // for (auto i = 0; i < src_fold; ++i) src.division_precision += src.division_precision;
    // auto ans = src.exp();
    // for (auto i = 0ull; i < cnt; ++i) ans *= ans;
    // ans.division_precision = src_prec;
    // ans.truncate();
    // cout << ans << endl;

    // auto fold = NEUNET_DEC_SEG_MAX;
    // auto cnt  = 0;
    // while (fold && !(fold % 2)) {
    //     fold /= 2;
    //     cout << ++cnt << ' ' << fold << ' ' << NEUNET_DEC_SEG_MAX - fold << endl;
    //     // std::printf("[%d][%llx][%llx]\n", ++cnt, fold, NEUNET_DEC_SEG_MAX - fold);
    // }
    // cout << endl;
    // for (auto i = 1; i <= cnt; ++i) {
    //     // std::printf("%llx\n", dec_bit_csh(i));
    //     cout << '[' << i << ']' << dec_bit_csh(i) << endl;
    // }

    // net_decimal test = 216,
    //             one  = 1;
    // one.division_precision = test.to_integer();
    // auto res  = one.exp(),
    //      base = res;
    // for (auto i = 1; i < test; ++i) res *= base;
    // res.truncate();
    // cout << res << endl;

    cout << '\n' << (NEUNET_CHRONO_TIME_POINT - ch_tm_pt) << "ms" << endl;
    return EXIT_SUCCESS;
}
