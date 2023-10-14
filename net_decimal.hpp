NEUNET_BEGIN

struct net_decimal_data final {
    /* 1|1010|0100|0131|0800|1920|0001.0001|0400|0021|7080|0400|005
     * {1,1920,800,131,100,1010,1}     {1,400,21,7080,400,50}
     */
    net_set<uint64_t> it, ft;

    void value_copy(const net_decimal_data &src) {
        it = src.it;
        ft = src.ft;
    }

    void value_move(net_decimal_data &&src) {
        it = std::move(src.it);
        ft = std::move(src.ft);
    }

    net_decimal_data() {}
    net_decimal_data(const net_decimal_data &src) { value_copy(src); }
    net_decimal_data(net_decimal_data &&src) { value_move(std::move(src)); }

    void reset() {
        it.reset();
        ft.reset();
    }

    ~net_decimal_data() { reset(); }

    net_decimal_data &operator=(const net_decimal_data &src) {
        value_copy(src);
        return *this;
    }
    net_decimal_data &operator=(net_decimal_data &&src) {
        value_move(std::move(src));
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const net_decimal_data &src) {
        if (src.it.length) for (auto i = src.it.length; i; --i) {
            auto pow_seg = NEUNET_DEC_SEG_MAX / 10,
                 seg_tmp = src.it[i - 1];
            if (i < src.it.length) while (pow_seg > seg_tmp) {
                pow_seg /= 10;
                os << 0;
            }
            os << seg_tmp;
        } else os << 0;
        if (!src.ft.length) return os;
        os << '.';
        for (auto i = 0ull; i < src.ft.length; ++i) {
            auto pow_seg = NEUNET_DEC_SEG_MAX / 10,
                 seg_tmp = src.ft[i];
            while (pow_seg > seg_tmp) {
                pow_seg /= 10;
                os << 0;
            }
            if ((i + 1) == src.ft.length) while (!(seg_tmp % 10)) seg_tmp /= 10;
            os << seg_tmp;
        }
        return os;
    }
};

uint64_t dec_dig_cnt(uint64_t seg, bool is_it) {
    auto ans = is_it ? 0ull : NEUNET_DEC_DIG_MAX;
    while (is_it ? seg : !(seg % 10)) {
        is_it ? ++ans : --ans;
        seg /= 10;
    }
    return ans;
}

uint64_t dec_dig_cnt(const net_decimal_data &src, bool is_it) {
    auto &seg = is_it ? src.it : src.ft;
    if (!seg.length) return 0;
    auto idx = seg.length - 1,
         ans = (idx) * NEUNET_DEC_DIG_MAX;
    return ans + dec_dig_cnt(seg[idx], is_it);
}

bool dec_is_zero(const net_decimal_data &src) { return !(src.it.length || src.ft.length); }

bool dec_is_one(const net_decimal_data &src) { return !src.ft.length && src.it.length == 1 && src.it[0] == 1; }

net_decimal_data dec_init(bool &sgn, long double src) {
    sgn = src < 0;
    src = std::abs(src);
    net_decimal_data ans;
    auto it_seg = (uint64_t)src;
    auto ft_seg = src - it_seg;
    if (it_seg > 0) {
        ans.it.init(1);
        ans.it[0] = it_seg;
    }
    if (!ft_seg) return ans;
    ans.ft.init(1);
    uint64_t dig_cnt = 0;
    while (dig_cnt++ < NEUNET_DEC_DIG_MAX) {
        if (ft_seg) ft_seg *= 10;
        ans.ft[0] *= 10;
        ans.ft[0] += ft_seg;
        if (ft_seg) ft_seg -= (int)ft_seg;
    }
    return ans;
}
net_decimal_data dec_init(bool &sgn, const std::string &src) {
    // verify
    auto str_len  = src.length(),
         dot_idx  = str_len,
         bgn_idx  = 0ull,
         end_idx  = dot_idx - 1;
    // symbol
    auto dot_flag = false;
    net_decimal_data ans;
    for (auto i = 0ull; i < str_len; ++i) if (src[i] < '0' || src[i] > '9') {
        if (src[i] == '.' && !dot_flag) {
            dot_flag = true;
            dot_idx  = i;
            continue;
        } 
        if (src[i] == '-' && !i) {
            sgn = true;
            ++bgn_idx;
            continue;
        }
        if (src[i] == '+' && !i) {
            ++bgn_idx;
            continue;
        }
        return ans;
    }
    // zero
    while (end_idx > dot_idx && src[end_idx] == '0') --end_idx;
    while (bgn_idx < dot_idx && src[bgn_idx] == '0') ++bgn_idx;
    if (end_idx == bgn_idx && dot_idx == bgn_idx) return ans;
    uint64_t seg_cnt = 0,
             seg_tmp = 0,
             tmp_cnt = 0;
    // integer
    if (bgn_idx < dot_idx) {
        tmp_cnt = dot_idx - bgn_idx;
        seg_cnt = tmp_cnt / NEUNET_DEC_DIG_MAX;
        if (tmp_cnt % NEUNET_DEC_DIG_MAX) ++seg_cnt;
        ans.it.init(seg_cnt, false);
        seg_cnt = 0;
        tmp_cnt = 1;
        for (auto i = dot_idx; i > bgn_idx; --i) {
            seg_tmp += (src[i - 1] - '0') * tmp_cnt;
            tmp_cnt *= 10;
            if (tmp_cnt == NEUNET_DEC_SEG_MAX) {
                ans.it[seg_cnt++] = seg_tmp;
                tmp_cnt           = 1;
                seg_tmp           = 0;
            }
        }
        if (seg_tmp) {
            ans.it[seg_cnt] = seg_tmp;
            seg_tmp     = 0;
        }
        seg_cnt = 0;
    }
    // float
    if (end_idx > dot_idx) {
        tmp_cnt   = end_idx - dot_idx;
        seg_cnt   = tmp_cnt / NEUNET_DEC_DIG_MAX;
        if (tmp_cnt % NEUNET_DEC_DIG_MAX) ++seg_cnt;
        ans.ft.init(seg_cnt, false);
        seg_cnt = 0;
        tmp_cnt = 0;
        for (auto i = dot_idx + 1; i <= end_idx; ++i) {
            seg_tmp *= 10;
            seg_tmp += src[i] - '0';
            ++tmp_cnt;
            if (tmp_cnt == NEUNET_DEC_DIG_MAX) {
                ans.ft[seg_cnt++] = seg_tmp;
                tmp_cnt           = 0;
                seg_tmp           = 0;
            }
        }
        if (seg_tmp) {
            while (tmp_cnt++ < NEUNET_DEC_DIG_MAX) seg_tmp *= 10;
            ans.ft[seg_cnt] = seg_tmp;
        }
    }
    return ans;
}

std::string dec_to_string(bool sgn, const net_decimal_data &src) {
    if (dec_is_zero(src)) return "0";
    std::string ans = sgn ? "-" : "";
    if (src.it.length) for (auto i = src.it.length; i; --i) {
        auto tmp = std::to_string(src.it[i - 1]);
        auto len = NEUNET_DEC_DIG_MAX - tmp.length();
        if (i < src.it.length) for (auto i = 0; i < len; ++i) ans.push_back('0');
        ans += tmp;
    } else ans += "0";
    if (src.ft.length) ans.push_back('.'); 
    for (auto i = 0ull; i < src.ft.length; ++i) {
        auto tmp = std::to_string(src.ft[i]);
        auto len = NEUNET_DEC_DIG_MAX - tmp.length();
        for (auto i = 0; i < len; ++i) ans.push_back('0');
        if ((i + 1) == src.ft.length) {
            len = tmp.length();
            while (tmp[len - 1] == '0') --len;
            tmp = tmp.substr(0, len);
        }
        ans += tmp;
    }
    return ans;
}
std::string dec_to_string(long long src) {
    std::string ans = "[";
    auto space = 4;
    if (src >= 0) --space;
    if (std::abs(src) < 100) --space;
    if (std::abs(src) < 10) --space;
    for (auto i = space; i < 4; ++i) ans.push_back(' ');
    ans += std::to_string(src);
    return ans + ']';
}

int64_t dec_to_integer(bool sgn, const net_decimal_data &src) {
    net_assert(src.it.length == 1, "neunet", "dec_to_integer(net_decimal_data)", "Value is out of limit.");
    if (dec_is_zero(src)) return 0;
    int64_t ans = src.it[0];
    return sgn ? -ans : ans;
}

long double dec_to_float(bool sgn, const net_decimal_data &src) {
    auto dig_cnt = dec_dig_cnt(src, true) + dec_dig_cnt(src, false);
    net_assert(dig_cnt <= NEUNET_DEC_VLD_DIG, "neunet", "dec_to_float(net_decimal_data)", "Value is greater than valid digit count.");
    if (dec_is_zero(src)) return 0;
    long double ans = 0;
    if (src.ft.length) {
        ans = src.ft[0];
        while (ans > 1) ans /= 10;
    }
    if (src.it.length) return ans + src.it[0];
    return ans;
}

net_decimal_data dec_float_part(net_decimal_data &integer_part, const net_decimal_data &src) {
    integer_part.it = src.it;
    integer_part.ft.reset();
    net_decimal_data ans;
    ans.ft = src.ft;
    return ans;
}

/* unsigned decimal digit part
 * NEUNET_DEC_CMP_EQL first is equal to second
 * NEUNET_DEC_CMP_LES first is less than second
 * NEUNET_DEC_CMP_GTR first is greater than second
 */
int dec_comp(const net_decimal_data &fst, const net_decimal_data &snd) {
    if (fst.it.length > snd.it.length) return NEUNET_DEC_CMP_GTR;
    if (fst.it.length < snd.it.length) return NEUNET_DEC_CMP_LES;
    for (auto i = fst.it.length; i; --i) {
        auto idx = i - 1;
        if (fst.it[idx] > snd.it[idx]) return NEUNET_DEC_CMP_GTR;
        if (fst.it[idx] < snd.it[idx]) return NEUNET_DEC_CMP_LES;
    }
    for (auto i = 0; i < fst.ft.length; ++i) {
        if (i == snd.ft.length) return NEUNET_DEC_CMP_GTR;
        if (fst.ft[i] > snd.ft[i]) return NEUNET_DEC_CMP_GTR;
        if (fst.ft[i] < snd.ft[i]) return NEUNET_DEC_CMP_LES;
    }
    if (fst.ft.length == snd.ft.length) return NEUNET_DEC_CMP_EQL;
    else return NEUNET_DEC_CMP_LES;
}
int dec_comp(bool sgn, int comp_res) {
    if (comp_res == NEUNET_DEC_CMP_GTR && sgn) return NEUNET_DEC_CMP_LES;
    if (comp_res == NEUNET_DEC_CMP_LES && sgn) return NEUNET_DEC_CMP_GTR;
    return comp_res;
}

uint64_t dec_sub(bool &carry, uint64_t minu, uint64_t subt) {
    auto cay = carry;
    carry    = subt > minu;
    if (carry) return NEUNET_DEC_SEG_MAX - subt + minu - cay;
    auto ans = minu - subt;
    if (cay) return dec_sub(carry, ans, 1);
    return ans;
}

// unsigned number digit segment
uint64_t dec_add(bool &carry, uint64_t fst, uint64_t snd) {
    auto dif = NEUNET_DEC_SEG_MAX - fst;
    auto cay = carry;
    carry    = dif <= snd;
    if (carry) return snd - dif + cay;
    auto ans = fst + snd;
    if (cay) return dec_add(carry, ans, 1);
    return ans;
}
uint64_t dec_add(bool &carry, uint64_t fst, uint64_t snd, bool subt) { return subt ? dec_sub(carry, fst, snd) : dec_add(carry, fst, snd); }
/* unsigned segment number
 * minuhend (first) segment should be greater than subtrahend (second) segment for true value of parameter subtract
 */
net_decimal_data dec_add(const net_decimal_data &fst, const net_decimal_data &snd, bool subt = false) {
    auto it_len  = std::max(fst.it.length, snd.it.length),
         ft_len  = std::max(fst.ft.length, snd.ft.length),
         ans_idx = it_len + ft_len;
    net_set<uint64_t> ans_tmp(ans_idx);
    net_decimal_data ans;
    auto carry = false;
    for (auto i = ft_len; i; --i) {
        auto idx_tmp = i - 1;
        auto seg_tmp = dec_add(carry,
                               fst.ft.length > idx_tmp ? fst.ft[idx_tmp] : 0,
                               snd.ft.length > idx_tmp ? snd.ft[idx_tmp] : 0,
                               subt);
        if (seg_tmp || ans_idx < ans_tmp.length) ans_tmp[--ans_idx] = seg_tmp;
    }
    if (ans_idx < ans_tmp.length) ans.ft = ans_tmp.sub_set(ans_idx, ans_tmp.length - 1);
    ans_idx = 0;
    for (auto i = 0ull; i < it_len; ++i) ans_tmp[ans_idx++] = dec_add(carry,
                                                                      fst.it.length > i ? fst.it[i] : 0,
                                                                      snd.it.length > i ? snd.it[i] : 0,
                                                                      subt);
    if (carry) ans_tmp[ans_idx++] = carry;
    else while (ans_idx && !ans_tmp[ans_idx - 1]) --ans_idx;
    if (ans_idx) ans.it = ans_tmp.sub_set(0, ans_idx - 1);
    return ans;
}

void dec_coe(uint64_t coe[]) {
    coe[0] %= NEUNET_DEC_MUL_POW;
    coe[1] %= NEUNET_DEC_MUL_SQR;
    coe[1] /= NEUNET_DEC_MUL_POW;
    coe[2] /= NEUNET_DEC_MUL_SQR;
}

uint64_t dec_mul_val(const net_decimal_data &src, uint64_t idx) {
    if (idx < src.ft.length) return src.ft[src.ft.length - idx - 1];
    else return src.it[idx - src.ft.length];
}

void dec_mul_carry(uint64_t &carry, bool &ca_add) { if (ca_add) {
    ca_add = false;
    ++carry;
} }

uint64_t dec_mul(uint64_t &carry, uint64_t fst, uint64_t snd) {
    /*
    9999999|999999|999999 ^ 2 =
    9999999999999999998|0000000000000000001
                               999998000001
                        1999996000002
                2099997|6000003
          1999997800000|2
    99999980000001
    */
    uint64_t fst_coe[3] = {fst, fst, fst},
             snd_coe[3] = {snd, snd, snd},
             ans_coe[5] = {0};
    dec_coe(fst_coe);
    dec_coe(snd_coe);
    carry = 0;
    for (auto i = 0; i < 3; ++i) for (auto j = 0; j < 3; ++j) ans_coe[i + j] += fst_coe[i] * snd_coe[j];
    auto ca_add = false;
    auto ans_1  = ans_coe[1] * NEUNET_DEC_MUL_POW,
         ans_2  = ans_coe[2] % NEUNET_DEC_MUL_END * NEUNET_DEC_MUL_SQR,
         ans_3  = ans_coe[3] % 10 * NEUNET_DEC_MUL_CUB,
         ans    = dec_add(ca_add, ans_coe[0], ans_1);
    dec_mul_carry(carry, ca_add);
    ans = dec_add(ca_add, ans, ans_2);
    dec_mul_carry(carry, ca_add);
    ans = dec_add(ca_add, ans, ans_3);
    dec_mul_carry(carry, ca_add);
    ans_coe[2] /= NEUNET_DEC_MUL_END;
    ans_coe[3] /= 10;
    carry      += ans_coe[2] + ans_coe[3] + ans_coe[4] * (NEUNET_DEC_MUL_POW / 10);
    return ans;
}
void dec_mul(net_set<uint64_t> &ans_coe, uint64_t fst, uint64_t snd, uint64_t fst_idx, uint64_t snd_idx) {
    uint64_t coe_idx = fst_idx + snd_idx,
             carry   = 0;
    auto     ca_add  = false;
    ans_coe[coe_idx] = dec_add(ca_add, ans_coe[coe_idx], dec_mul(carry, fst, snd));
    carry           += ca_add;
    ca_add           = false;
    ++coe_idx;
    ans_coe[coe_idx] = dec_add(ca_add, ans_coe[coe_idx], carry);
    while (ca_add) {
        ++coe_idx;
        ans_coe[coe_idx] = dec_add(ca_add, ans_coe[coe_idx], 0);
    }
}
net_decimal_data dec_mul(const net_decimal_data &fst, const net_decimal_data &snd) {
    auto fst_len = fst.it.length + fst.ft.length,
         snd_len = snd.it.length + snd.ft.length,
         ft_len  = fst.ft.length + snd.ft.length;
    net_set<uint64_t> ans_coe(fst_len * snd_len * 2);
    for (auto i = 0ull; i < fst_len; ++i) {
        auto fst_val = dec_mul_val(fst, i);
        for (auto j = 0ull; j < snd_len; ++j) {
            auto snd_val = dec_mul_val(snd, j);
            dec_mul(ans_coe, fst_val, snd_val, i, j);
        }
    }
    net_decimal_data ans;
    auto idx = 0ull;
    while (idx < ft_len && !ans_coe[idx]) ++idx;
    if (idx < ft_len) {
        ans.ft = ans_coe.sub_set(idx, ft_len - 1);
        ans.ft.reverse();
    }
    idx = ans_coe.length;
    while (ft_len < idx && !ans_coe[idx - 1]) --idx;
    if (ft_len < idx) ans.it = ans_coe.sub_set(ft_len, idx - 1);
    return ans;
}

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
    else if (ft_cnt) dig_cnt = NEUNET_DEC_DIG_MAX;
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

uint8_t dec_div(net_set<uint8_t> &divd, const net_set<uint8_t> &divr) {
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
    if (ans_coe >= 10) ans_coe = 9;
    while (true) {
        auto carry = 0ll;
        
        // for (auto i = 0ull; i < divd.length; ++i) std::cout << dec_to_string(divd[i]); std::cout << '\n';
        // for (auto i = 0ull; i < divr.length; ++i) std::cout << dec_to_string(ans_coe * divr[i]); std::cout << '\n';
        
        if (divd.length > divr.length) divd_tmp.copy(divr.length - 1, divd, divr.length, divd.length - divr.length);
        for (auto i = divr.length; i > 1; --i) {
            auto idx = i - 1;
            if (idx < divd.length) divd_tmp[idx - 1] = dec_carry(carry, divd[idx] - ans_coe * divr[idx]);
            else divd_tmp[idx - 1] = dec_carry(carry, 0 - divr[idx]);
        }
        carry += divd[0] - ans_coe * divr[0];

        // std::cout << dec_to_string(carry);
        // for (auto i = 0ull; i < divd_tmp.length; ++i) std::cout << dec_to_string(divd_tmp[i]); std::cout << '\n';

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
net_decimal_data dec_div(const net_decimal_data &divd, const net_decimal_data &divr, uint64_t prec) {
    net_assert(!dec_is_zero(divr), "neunet", "dec_div(net_decimal_data)", "Divisor could not be 0.");
    auto divd_ft_cnt = dec_dig_cnt(divd, false),
         divd_it_cnt = dec_dig_cnt(divd, true),
         divr_ft_cnt = dec_dig_cnt(divr, false),
         divr_it_cnt = dec_dig_cnt(divr, true),
         ans_it_cnt  = 1ull,
         ans_idx_tmp = 0ull;
    net_set<uint8_t> divr_coe(divr_it_cnt + divr_ft_cnt);
    net_set<uint8_t> divd_coe(divd_it_cnt + divd_ft_cnt);
    dec_coe(divd_coe, divd);
    dec_coe(divr_coe, divr);
    if (divd_it_cnt < divr_it_cnt) ans_idx_tmp = divr_it_cnt - divd_it_cnt;
    else ans_it_cnt = divd_it_cnt - divr_it_cnt + 1;
    net_set<uint8_t> ans_coe(ans_it_cnt + prec + 2);
    for (auto i = 0; i < ans_coe.length; ++i) ans_coe[i] = dec_div(divd_coe, divr_coe);

    // std::cout << "ans" << std::endl;
    // for (int tmp : ans_coe) std::cout << tmp << std::endl;

    return dec_coe(ans_coe, prec + 2);
}

// return quotient
net_decimal_data dec_rem(net_decimal_data &divd_rem, const net_decimal_data &divr) {
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

bool dec_gcd(net_decimal_data &fst, net_decimal_data &snd) {
    auto comp = dec_comp(fst, snd);
    if (comp == NEUNET_DEC_CMP_EQL) return false;
    if (comp == NEUNET_DEC_CMP_GTR) {
        dec_rem(fst, snd);
        // std::cout << dec_to_string(false, fst) << '\n' << dec_to_string(false, snd) << '\n' << std::endl;
        if (dec_is_zero(fst)) return true;
    }
    if (comp == NEUNET_DEC_CMP_LES) {
        dec_rem(snd, fst);
        // std::cout << dec_to_string(false, fst) << '\n' << dec_to_string(false, snd) << '\n' << std::endl;
        if (dec_is_zero(snd)) return false;
    }
    return dec_gcd(fst, snd);
}

net_decimal_data dec_e10(net_decimal_data &e10, const net_decimal_data &src) {
    auto sgn = false;
    if (!src.ft.length) {
        e10 = dec_init(sgn, 1);
        return src;
    }
    auto it_seg = src.ft;
    it_seg.reverse();
    it_seg = it_seg.unit(src.it);
    /*
    41534 3607414684000115610.1056011112001114511 000564861
    41534|3607414684000115610|1056011112001114511|   564861000000000
    */
    net_decimal_data ans;
    auto seg_tmp = it_seg[0],
         pow_rem = 1ull,
         e10_seg = src.ft.length;
    if (seg_tmp % 10) {
        // mice segment for the multiple of 19
        e10.it.init(e10_seg + 1);
        e10.it[e10_seg] = 1;
        ans.it = std::move(it_seg);
        return ans;
    } else e10.it.init(e10_seg--);
    while (!(seg_tmp % 10)) {
        // 564861000000000 -> 564861000000000 / 1000000000
        pow_rem *= 10;
        seg_tmp /= 10;
    }
    auto pow_cnt = NEUNET_DEC_SEG_MAX / pow_rem,
         end_idx = it_seg.length - 1;
    /*
    10000000000000000000
              1000000000
     1000000000
    */
    for (auto i = 0ull; i < end_idx; ++i) {
        /*
        1056011112001114511|   564861000000000
        564861000000000 / 1000000000 = 564861
        1056011112001114511 % 1000000000 * 1000000000 = 2001114511000000000
        2001114511000000000 + 564861 = 2001114511000564861
        1056011112001114511 / 1000000000 = 105601111
        105601111 | 2001114511000564861
        */
        it_seg[i] /= pow_rem;
        it_seg[i] += it_seg[i + 1] % pow_rem * pow_cnt;
    }
    e10.it[e10_seg]  = pow_cnt;
    it_seg[end_idx] /= pow_rem;
    while (!it_seg[end_idx]) --end_idx;
    ans.it = it_seg.sub_set(0, end_idx);
    return ans;
}

void dec_truncate(net_decimal_data &src, uint64_t prec) {
    auto dig_cnt = dec_dig_cnt(src, false);
    auto sgn_tmp = false;
    if (prec >= dig_cnt) return;
    net_decimal_data carry;
    if (prec) {
        auto carry_seg = NEUNET_DEC_SEG_MAX / 10,
             seg_cnt   = 1ull;
        for (auto i = 0ull; i < prec; ++i) if (carry_seg == 1) {
            carry_seg = NEUNET_DEC_SEG_MAX / 10;
            ++seg_cnt;
        } else carry_seg /= 10;
        carry.ft.init(seg_cnt, false);
        carry.ft[seg_cnt - 1] = carry_seg;
    } else carry = dec_init(sgn_tmp, 1);
    prec += 2;
    uint64_t tgt_idx = prec / NEUNET_DEC_DIG_MAX,
             seg_len = prec % NEUNET_DEC_DIG_MAX,
             dig_idx = NEUNET_DEC_DIG_MAX;
    if (seg_len) ++tgt_idx;
    else seg_len = dig_idx;
    src.ft.init(tgt_idx--);
    auto curr_seg = src.ft[tgt_idx],
         seg_pow  = 1ull,
         last_dig = 0ull;
    // get last 2 digits after accuracy digits
    while (dig_idx-- > seg_len) {
        last_dig             = seg_pow * (curr_seg % 10);
        curr_seg            /= 10;
        seg_pow             *= 10;
        src.ft[tgt_idx] -= last_dig;
    }
    last_dig = curr_seg % 10;
    curr_seg /= 10;
    // last digit of the 2;
    if (last_dig >= 5) {
        // carry
        src = dec_add(src, carry);
        ++curr_seg;
    }
    last_dig *= seg_pow;
    seg_pow  *= 10;
    if (dig_idx) {
        // not the end of current segment
        src.ft[tgt_idx] -= last_dig;
        last_dig         = curr_seg % 10;
        if (last_dig < 5) src.ft[tgt_idx] -= last_dig * seg_pow;
        if (!src.ft[tgt_idx]) src.ft.init(tgt_idx);
    } else {
        // current segment ends
        src.ft.init(tgt_idx--);
        last_dig = src.ft[tgt_idx] % 10;
        if (last_dig < 5) src.ft[tgt_idx] -= last_dig;
    }
}

net_decimal_data dec_add(bool &ans_sgn, const net_decimal_data &fst, bool fst_sgn, const net_decimal_data &snd, bool snd_sgn) {
    if (dec_is_zero(fst)) {
        ans_sgn = snd_sgn;
        return snd;
    }
    if (dec_is_zero(snd)) {
        ans_sgn = fst_sgn;
        return fst;
    }
    ans_sgn = false;
    if (fst_sgn != snd_sgn) {
        auto cmp_res = dec_comp(snd, fst);
        if (cmp_res == NEUNET_DEC_CMP_GTR) {
            if (snd_sgn) ans_sgn = true;
            return dec_add(snd, fst, true);
        }
        if (cmp_res == NEUNET_DEC_CMP_LES) {
            if (fst_sgn) ans_sgn = true;
            return dec_add(fst, snd, true);
        }
        return {};
    }
    if (fst_sgn) ans_sgn = true;
    return dec_add(fst, snd);
}

net_decimal_data dec_sub(bool &ans_sgn, const net_decimal_data &minu, bool minu_sgn, const net_decimal_data &subt, bool subt_sgn) {
    if (dec_is_zero(subt)) {
        ans_sgn = minu_sgn;
        return minu;
    }
    if (dec_is_zero(minu)) {
        ans_sgn = !minu_sgn;
        return subt;
    }
    ans_sgn = false;
    if (minu_sgn == subt_sgn) {
        auto cmp_res = dec_comp(minu, subt);
        if (cmp_res == NEUNET_DEC_CMP_GTR) {
            if (minu_sgn) ans_sgn = true;
            return dec_add(minu, subt, true);
        }
        if (cmp_res == NEUNET_DEC_CMP_LES) {
            if (!minu_sgn) ans_sgn = true;
            return dec_add(subt, minu, true);
        }
        return {};
    }
    if (minu_sgn) ans_sgn = true;
    return dec_add(minu, subt);
}

net_decimal_data dec_mul(bool &ans_sgn, const net_decimal_data &fst, bool fst_sgn, const net_decimal_data &snd, bool snd_sgn) {
    ans_sgn = fst_sgn != snd_sgn;
    if (dec_is_zero(fst) || dec_is_zero(snd)) return {};
    if (dec_is_one(fst)) return snd;
    if (dec_is_one(snd)) return fst;
    return dec_mul(fst, snd);
}

net_decimal_data dec_div(bool &ans_sgn, const net_decimal_data &divd, bool divd_sgn, const net_decimal_data &divr, bool divr_sgn, uint64_t prec) {
    ans_sgn = divd_sgn != divr_sgn;
    if (dec_is_zero(divd)) return {};
    if (dec_is_zero(divr)) return divd;
    return dec_div(divd, divr, prec);
}

struct net_decimal_frac final {
    bool red = false;
    // numerator / denominator
    net_decimal_data num, den;

    void value_assign(const net_decimal_frac &src) { red = src.red; }

    void value_copy(const net_decimal_frac &src) {
        value_assign(src);
        num = src.num;
        den = src.den;
    }

    void value_move(net_decimal_frac &&src) {
        value_assign(src);
        num = std::move(src.num);
        den = std::move(src.den);
    }

    net_decimal_frac() {}
    net_decimal_frac(const net_decimal_frac &src) { value_copy(src); }
    net_decimal_frac(net_decimal_frac &&src) { value_move(std::move(src)); }

    void reset() {
        num.reset();
        den.reset();
        red = false;
    }

    ~net_decimal_frac() { reset(); }

    net_decimal_frac &operator=(const net_decimal_frac &src) {
        value_copy(src);
        return *this;
    }
    net_decimal_frac &operator=(net_decimal_frac &&src) {
        value_move(std::move(src));
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const net_decimal_frac &src) {
        os << "[Numerator][\n";
        os << src.num << '\n';
        os << "][Denominator][\n";
        os << src.den << '\n]';
    }
};

bool dec_is_zero(net_decimal_frac &src) {
    auto ans = dec_is_zero(src.num);
    if (!ans) return false;
    src.den = dec_init(src.red, 1);
    src.red = true;
    return true;
}

bool dec_is_one(net_decimal_frac &src) {
    if (dec_comp(src.den, src.num) != NEUNET_DEC_CMP_EQL) return false;
    src.num = dec_init(src.red, 1);
    src.den = dec_init(src.red, 1);
    src.red = true;
    return true;
}

bool dec_is_valid(const net_decimal_frac &src) { return !dec_is_zero(src.den); }

void dec_norm(net_decimal_frac &src) {
    net_decimal_data e10;
    if (src.num.ft.length) {
        src.num = dec_e10(e10, src.num);
        src.den = dec_mul(e10, src.den);
    }
    if (src.den.ft.length) {
        src.den = dec_e10(e10, src.den);
        src.num = dec_mul(e10, src.num);
    }
}

// return compare result of the numerator and denominator
int dec_reduct(net_decimal_frac &src) {
    if (src.red) return -1;
    else src.red = true;
    if (dec_is_one(src)) return NEUNET_DEC_CMP_EQL;
    if (dec_is_zero(src)) return NEUNET_DEC_CMP_LES;
    auto ans = dec_comp(src.num, src.den);
    if (ans == NEUNET_DEC_CMP_EQL) {
        auto sgn = false;
        src.num  = dec_init(sgn, 1);
        src.den  = dec_init(sgn, 1);
        return ans;
    }
    dec_norm(src);
    auto fst = src.num,
         snd = src.den;
    auto gcd_id = dec_gcd(fst, snd);
    net_decimal_data gcd_val;
    if (gcd_id) gcd_val = std::move(snd);
    else gcd_val = std::move(fst);
    src.num = dec_rem(src.num, gcd_val);
    src.den = dec_rem(src.den, gcd_val);
    return ans;
}

net_decimal_frac dec_init(const net_decimal_data &src, bool reduct = false) {
    net_decimal_frac ans;
    ans.den = dec_init(ans.red, 1);
    ans.num = src;
    if (reduct) dec_reduct(ans);
    return ans;
}
net_decimal_frac dec_init(const net_decimal_data &num, const net_decimal_data &den, bool reduct = false) {
    net_assert(!dec_is_zero(den), "neunet", "dec_init(net_decimal_data)", "Denominator should not be 0.");
    net_decimal_frac ans;
    ans.num = num;
    ans.den = den;
    if (reduct) dec_reduct(ans);
    return ans;
}

std::string dec_to_string(bool sgn, const net_decimal_frac &src, uint64_t prec) {
    if (dec_is_one(src.den)) return dec_to_string(sgn, src.num);
    return dec_to_string(sgn, dec_div(src.num, src.den, prec));
}

int64_t dec_to_integer(bool sgn, net_decimal_frac &src) {
    dec_reduct(src);
    if (dec_is_one(src.den)) return dec_to_integer(sgn, src.num);
    net_decimal_data num = src.num;
    num = dec_rem(num, src.den);
    return dec_to_integer(sgn, num);
}

long double dec_to_float(bool sgn, net_decimal_frac &src) {
    dec_reduct(src);
    if (dec_is_one(src.den)) return dec_to_float(sgn, src.num);
    auto ans = dec_div(src.num, src.den, NEUNET_DEC_VLD_DIG);
    auto len = dec_dig_cnt(ans, true);
    net_assert(len <= NEUNET_DEC_VLD_DIG, "neunet", "dec_to_float(net_decimal_frac)", "Value is greater than valid digit count.");
    if (len == NEUNET_DEC_VLD_DIG) ans.ft.reset();
    else dec_truncate(ans, NEUNET_DEC_VLD_DIG - len - 2);
    return dec_to_float(sgn, ans);
}

net_decimal_frac dec_float_part(net_decimal_data &integer_part, net_decimal_frac &src) {
    auto cmp_res = dec_reduct(src);
    net_decimal_frac ans;
    if (dec_is_one(src.den)) {
        integer_part  = src.num;
        ans.den       = dec_init(ans.red, 1);
        ans.red       = true;
        return ans;
    }
    ans = src;
    if (cmp_res == NEUNET_DEC_CMP_LES) {
        integer_part.reset();
        return ans;
    }
    integer_part = dec_rem(ans.num, ans.den);
    return ans;
}

net_decimal_frac dec_add(bool &ans_sgn, net_decimal_frac &fst, bool fst_sgn, net_decimal_frac &snd, bool snd_sgn, bool reduct) {
    net_decimal_frac ans;
    if (dec_comp(fst.den, snd.den) == NEUNET_DEC_CMP_EQL) {
        ans.num = dec_add(ans_sgn, fst.num, fst_sgn, snd.num, snd_sgn);
        ans.den = fst.den;
    } else {
        if (reduct) {
            dec_reduct(fst);
            dec_reduct(snd);
        }
        ans.den = dec_mul(fst.den, snd.den);
        ans.num = dec_add(ans_sgn, dec_mul(fst.num, snd.den), fst_sgn, dec_mul(snd.num, fst.den), snd_sgn);
    }
    if (reduct) dec_reduct(ans);
    return ans;
}

net_decimal_frac dec_sub(bool &ans_sgn, net_decimal_frac &minu, bool minu_sgn, net_decimal_frac &subt, bool subt_sgn, bool reduct) {
    net_decimal_frac ans;
    if (dec_comp(minu.den, subt.den) == NEUNET_DEC_CMP_EQL) {
        ans.num = dec_sub(ans_sgn, minu.num, minu_sgn, subt.num, subt_sgn);
        ans.den = minu.den;
    } else {
        if (reduct) {
            dec_reduct(minu);
            dec_reduct(subt);
        }
        ans.den = dec_mul(minu.den, subt.den);
        ans.num = dec_sub(ans_sgn, dec_mul(minu.num, subt.den), minu_sgn, dec_mul(subt.num, minu.den), subt_sgn);
    }
    if (reduct) dec_reduct(ans);
    return ans;
}

net_decimal_frac dec_mul(bool &ans_sgn, net_decimal_frac &fst, bool fst_sgn, net_decimal_frac &snd, bool snd_sgn, bool reduct) {
    if (reduct) {
        dec_reduct(fst);
        dec_reduct(snd);
    }
    net_decimal_frac ans;
    ans.num = dec_mul(ans_sgn, fst.num, true, snd.num, true);
    ans.den = dec_mul(ans_sgn, fst.den, true, snd.den, true);
    ans_sgn = fst_sgn != snd_sgn;
    if (reduct) dec_reduct(ans);
    return ans;
}

net_decimal_frac dec_div(bool &ans_sgn, net_decimal_frac &divd, bool divd_sgn, net_decimal_frac &divr, bool divr_sgn, bool reduct) {
    if (reduct) {
        dec_reduct(divd);
        dec_reduct(divr);
    }
    net_decimal_frac ans;
    ans.num = dec_mul(ans_sgn, divd.num, true, divr.den, true);
    ans.den = dec_mul(ans_sgn, divd.den, true, divr.num, true);
    ans_sgn = divd_sgn != divr_sgn;
    if (reduct) dec_reduct(ans);
    return ans;
}

class net_decimal {
public:
    bool modulus = false;

    static std::atomic_uint64_t division_precision;

    static std::atomic_bool fraction_reduct;

protected:
    void value_assign(const net_decimal &src) {
        sgn     = src.sgn;
        dec     = src.dec;
        modulus = src.modulus;
    }

    void value_copy(const net_decimal &src) {
        value_assign(src);
        base = src.base;
        frac = src.frac;
    }

    void value_move(net_decimal &&src) {
        value_assign(src);
        base = std::move(src.base);
        frac = std::move(src.frac);
    }

    bool is_zero() {
        if (dec) return dec_is_zero(base);
        if (!dec_is_zero(frac)) return false;
        base = dec_init(dec, 0);
        dec  = true;
        return true;
    }

    bool is_one() {
        if (dec) return dec_is_one(base);
        if (!dec_is_one(frac)) return false;
        base = dec_init(dec, 1);
        dec  = true;
        return true;
    }

    static net_decimal ln_4(bool show_iter = false) {
        net_decimal b {1},
                    c {2},
                    s {"0.6"},
                    o {"0.36"},
                    p {0},
                    q {b};
        net_decimal_data m, n;
        do {
            if (show_iter) std::cout << n << std::endl;
            m = std::move(n);
            if (b == q) p += s;
            else {
                p  = b * p + s * q;
                q *= b;
            }
            s *= o;
            b += c;
            n  = dec_div(p.base, q.base, division_precision);
        } while (dec_comp(m, n) != NEUNET_DEC_CMP_EQL);
        net_decimal ans;
        ans.sgn  = p.sgn != q.sgn;
        ans.base = std::move(n);
        return 2 * ans;
    }
    // ln(x + 1) = (-1 < x <= 1)
    net_decimal ln_2(bool show_iter = false) {
        net_decimal b {1},
                    p {0},
                    q {b},
                    x {*this - 1},
                    o {x},
                    c {b};
        net_decimal_data m, n;
        do {
            if (show_iter) std::cout << n << std::endl;
            m = std::move(n);
            if (b == q) p += x;
            else {
                p  = b * p + c * x * q;
                q *= b;
            }
            x *= o;
            n  = dec_div(p.base, q.base, division_precision);
            ++b;
            c.sgn = !c.sgn;
        } while (dec_comp(m, n) != NEUNET_DEC_CMP_EQL);
        net_decimal ans;
        ans.base = std::move(n);
        ans.sgn  = p.sgn != q.sgn;
        return ans;
    }
    net_decimal ln_proto(bool show_iter = false) {
        net_decimal b {1},
                    p {0},
                    q {b},
                    o {p},
                    u {*this - 1},
                    v {*this + 1},
                    h {u * u},
                    s {v * v};
        net_decimal_data m, n;
        do {
            if (show_iter) std::cout << n << std::endl;
            m = std::move(n);
            if (q == v) p += u;
            else {
                o  = b * v;
                p  = o * p + u * q;
                q *= o;
            }
            u *= h;
            v *= s;
            b += 2;
            n  = dec_div(p.base, q.base, division_precision);
        } while (dec_comp(m, n) != NEUNET_DEC_CMP_EQL);
        net_decimal ans;
        ans.base = std::move(n);
        ans.sgn  = p.sgn != q.sgn;
        return 2 * ans;
    }

public:
    net_decimal() {}
    net_decimal(const net_decimal &src) { value_copy(src); }
    net_decimal(net_decimal &&src) { value_move(std::move(src)); }
    callback_arg net_decimal(const arg &src) {
        static_assert(neunet_dec_init_expr);
        base = dec_init(sgn, src);
    }

    std::string to_string() const {
        if (dec) return dec_to_string(sgn, base);
        return dec_to_string(sgn, frac, division_precision);
    }

    int64_t to_integer() {
        if (dec) return dec_to_integer(sgn, base);
        return dec_to_integer(sgn, frac);
    }

    long double to_float() {
        if (dec) return dec_to_float(sgn, base);
        return dec_to_float(sgn, frac);
    }

    net_decimal float_part(net_decimal &integer_part) {
        integer_part.reset();
        integer_part.sgn = sgn;
        net_decimal ans;
        ans.sgn = sgn;
        if (dec) ans.base = dec_float_part(integer_part.base, base);
        else {
            ans.frac = dec_float_part(integer_part.base, frac);
            ans.dec  = false;
        }
        return ans;
    }    

    int comp(net_decimal &src) {
        if (sgn != src.sgn) return sgn ? NEUNET_DEC_CMP_LES : NEUNET_DEC_CMP_GTR;
        if (dec && dec == src.dec) return dec_comp(sgn, dec_comp(base, src.base));
        if (!dec_is_valid(frac)) frac = dec_init(base);
        if (!dec_is_valid(src.frac)) src.frac = dec_init(src.base);
        auto fst_num = dec_mul(frac.num, src.frac.den),
             snd_num = dec_mul(src.frac.num, frac.den);
        return dec_comp(fst_num, snd_num);
    }

    void reset() {
        base.reset();
        frac.reset();
        dec = true;
        sgn = false;
    }

    ~net_decimal() { reset(); }

    static net_decimal pi() {
        net_decimal p {0},
                    q {1},
                    a {0.25}, // +2
                    b {2},    // +4
                    c {5},    // +8
                    d {6},    // +8
                    o {0.0625},
                    s {q};
        for (auto i = 0ull; i < division_precision; ++i) {
            net_decimal u = 1,
                        v = a;
            u  = b * u - v;
            v *= b;
            u  = c * u - v;
            v *= c;
            u  = d * u - v;
            v *= d;
            u *= s;
            if (q == v) p += u;
            else {
                p  = v * p + u * q;
                q *= v;
            }
            s *= o;
            a += 2;
            b += 4;
            c += 8;
            d += 8;
        }
        return p / q;
    }

    net_decimal ln() {
        if (is_zero() || is_one() || sgn) return 0.;
        net_decimal cnt  {0},
                    base {*this};
        while (base > 1) {
            base *= 0.25;
            ++cnt;
        }
        return base.ln_2() + cnt * ln_4();
    }

    net_decimal exp() {
        net_decimal c {1},
                    p {c},
                    q {c},
                    v {c},
                    u {*this};
        net_decimal_data m, n;
        do {
            m = std::move(n);
            if (q == v) p += u;
            else {
                p  = v * p + u * q;
                q *= v;
            }
            ++c;
            u *= *this;
            v *= c;
            n  = dec_div(p.base, q.base, division_precision);
        } while (dec_comp(m, n) != NEUNET_DEC_CMP_EQL);
        net_decimal ans;
        ans.base = std::move(n);
        ans.sgn  = p.sgn != q.sgn;
        return ans;
    }
    
protected:
    net_decimal_data base;

    net_decimal_frac frac;

    bool dec = true,
         sgn = false;

public:
    net_decimal &operator=(const net_decimal &src) {
        value_copy(src);
        return *this;
    }
    net_decimal &operator=(net_decimal &&src) {
        value_move(std::move(src));
        return *this;
    }

    callback_dec_s friend auto operator<=>(fst_dec_t &&fst, snd_dec_t &&snd) {
        neunet_type_if (neunet_dec_num(fst_dec_t)) return net_decimal {fst} <=> snd;
        neunet_type_elif (neunet_dec_num(snd_dec_t)) return fst <=> net_decimal {snd};
        neunet_type_else
        if (fst.sgn != snd.sgn) return fst.sgn ? std::strong_ordering::less : std::strong_ordering::greater;
        auto cmp_res = fst.comp(snd);
        if (cmp_res == NEUNET_DEC_CMP_LES) return fst.sgn ? std::strong_ordering::greater : std::strong_ordering::less;
        if (cmp_res == NEUNET_DEC_CMP_GTR) return fst.sgn ? std::strong_ordering::less : std::strong_ordering::greater;
        return std::strong_ordering::equal;
        neunet_type_endif
    }
    callback_dec_s friend bool operator==(fst_dec_t &&fst, snd_dec_t &&snd) {
        neunet_type_if (neunet_dec_num(fst_dec_t)) return net_decimal {fst} == snd;
        neunet_type_elif (neunet_dec_num(snd_dec_t)) return fst == net_decimal {snd};
        neunet_type_else return fst.comp(snd) == NEUNET_DEC_CMP_EQL;
        neunet_type_endif
    }

    net_decimal operator+() { return *this; }
    callback_dec_s friend net_decimal operator+(fst_dec_t &&fst, snd_dec_t &&snd) {
        // first is arithmetic type
        neunet_type_if (neunet_dec_num(fst_dec_t)) return net_decimal {fst} + snd;
        // second is arithmetic type
        neunet_type_elif (neunet_dec_num(snd_dec_t)) return fst + net_decimal {snd};
        // both are decimal
        neunet_type_else
        net_decimal ans;
        // TODO define addition function
        if (fst.dec && snd.dec) {
            ans.base = dec_add(ans.sgn, fst.base, fst.sgn, snd.base, snd.sgn);
            return ans;
        }
        if (!dec_is_valid(fst.frac)) fst.frac = dec_init(fst.base, fraction_reduct);
        if (!dec_is_valid(snd.frac)) snd.frac = dec_init(snd.base, fraction_reduct);
        ans.dec  = false;
        ans.frac = dec_add(ans.sgn, fst.frac, fst.sgn, snd.frac, snd.sgn, fraction_reduct);
        return ans;
        neunet_type_endif
    }
    callback_dec void operator+=(dec_t &&src) {
        neunet_type_if (neunet_dec_num(dec_t)) *this = *this + net_decimal {src};
        neunet_type_else *this = *this + src;
        neunet_type_endif
    }
    callback_dec_arg_ref friend void operator+=(arg &fst, net_decimal &snd) { fst += snd.to_float(); }
    net_decimal &operator++() {
        *this += 1;
        return *this;
    }
    net_decimal operator++(int) {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    net_decimal operator-() { return 0 - *this; }
    callback_dec_s friend net_decimal operator-(fst_dec_t &&minu, snd_dec_t &&subt) {
        // first is arithmetic type
        neunet_type_if (neunet_dec_num(fst_dec_t)) return net_decimal {minu} - subt;
        // second is arithmetic type
        neunet_type_elif (neunet_dec_num(snd_dec_t)) return minu - net_decimal {subt};
        // both are decimal
        neunet_type_else
        net_decimal ans;
        // TODO define subtraction function
        if (minu.dec && subt.dec) {
            ans.base = dec_sub(ans.sgn, minu.base, minu.sgn, subt.base, subt.sgn);
            return ans;
        }
        if (!dec_is_valid(minu.frac)) minu.frac = dec_init(minu.base, fraction_reduct);
        if (!dec_is_valid(subt.frac)) subt.frac = dec_init(subt.base, fraction_reduct);
        ans.dec  = false;
        ans.frac = dec_sub(ans.sgn, minu.frac, minu.sgn, subt.frac, subt.sgn, fraction_reduct);
        return ans;
        return ans;
        neunet_type_endif
    }
    callback_dec void operator-=(dec_t &&src) {
        neunet_type_if (neunet_dec_num(dec_t)) *this = *this - net_decimal {src};
        neunet_type_else *this = *this - src;
        neunet_type_endif
    }
    callback_dec_arg_ref friend void operator-=(arg &fst, net_decimal &snd) { fst -= snd.to_float(); }
    net_decimal &operator--() {
        *this -= 1;
        return *this;
    }
    net_decimal operator--(int) {
        auto tmp = *this;
        --(*this);
        return tmp;
    }

    callback_dec_s friend net_decimal operator*(fst_dec_t &&fst, snd_dec_t &&snd) {
        // first is arithmetic type
        neunet_type_if (neunet_dec_num(fst_dec_t)) return net_decimal {fst} * snd;
        // second is arithmetic type
        neunet_type_elif (neunet_dec_num(snd_dec_t)) return fst * net_decimal {snd};
        // both are decimal
        neunet_type_else
        net_decimal ans;
        ans.sgn = fst.sgn != snd.sgn;
        if (fst.dec && snd.dec) {
            ans.base = dec_mul(fst.base, snd.base);
            return ans;
        }
        if (!dec_is_valid(fst.frac)) fst.frac = dec_init(fst.base, fraction_reduct);
        if (!dec_is_valid(snd.frac)) snd.frac = dec_init(snd.base, fraction_reduct);
        ans.dec  = false;
        ans.frac = dec_mul(ans.sgn, fst.frac, fst.sgn, snd.frac, snd.sgn, fraction_reduct);
        return ans;
        neunet_type_endif
    }
    callback_dec void operator*=(dec_t &&src) {
        neunet_type_if (neunet_dec_num(dec_t)) *this = *this * net_decimal {src};
        neunet_type_else *this = *this * src;
        neunet_type_endif
    }
    callback_dec_arg_ref friend void operator*=(arg &fst, net_decimal &snd) { fst *= snd.to_float(); }

    callback_dec_s friend net_decimal operator/(fst_dec_t &&divd, snd_dec_t &&divr) {
        // first is arithmetic type
        neunet_type_if (neunet_dec_num(fst_dec_t)) return net_decimal {divd} / divr;
        // second is arithmetic type
        neunet_type_elif (neunet_dec_num(snd_dec_t)) return divd / net_decimal {divr};
        // both are decimal
        neunet_type_else
        net_decimal ans;
        ans.sgn = divr.sgn != divd.sgn;
        if (divr.is_one()) {
            ans = divd;
            return ans;
        }
        ans.dec = false;
        if (divd.dec && divd.dec) {
            ans.frac = dec_init(divd.base, divr.base, fraction_reduct);
            return ans;
        }
        if (!dec_is_valid(divd.frac)) divd.frac = dec_init(divd.base, fraction_reduct);
        if (!dec_is_valid(divr.frac)) divr.frac = dec_init(divr.base, fraction_reduct);
        ans.frac = dec_div(ans.sgn, divd.frac, divd.sgn, divr.frac, divr.sgn, fraction_reduct);
        return ans;
        neunet_type_endif
    }
    callback_dec void operator/=(dec_t &&divr) {
        neunet_type_if (neunet_dec_num(dec_t)) *this = *this / net_decimal {divr};
        neunet_type_else *this = *this / divr;
        neunet_type_endif
    }
    callback_dec_arg_ref friend void operator/=(arg &fst, net_decimal &snd) { fst /= snd.to_float(); }

    callback_dec_s friend net_decimal operator%(fst_dec_t &&divd, snd_dec_t &&divr) {
        // first is arithmetic type
        neunet_type_if (neunet_dec_num(fst_dec_t)) return net_decimal {divd} % divr;
        // second is arithmetic type
        neunet_type_elif (neunet_dec_num(snd_dec_t)) return divd % net_decimal {divr};
        // both are decimal
        neunet_type_else
        net_decimal ans;
        ans.sgn = divd.sgn;
        if (divd.dec) ans.base = divd.base;
        else {
            auto tmp = dec_float_part(ans.base, divd.frac);
            net_assert(!dec_is_zero(tmp), "net_decimal", "%", "Dividend could not be 0.");
        }
        if (divr.dec) {
            dec_rem(ans.base, divr.base);
            if (divr.modulus && divd.sgn != divr.sgn) ans.base = dec_add(ans.sgn, ans.base, false, divr.base, divr.sgn);
        } else {
            net_decimal_data divr_base;
            auto tmp = dec_float_part(divr_base, divr.frac);
            net_assert(!dec_is_zero(tmp), "net_decimal", "%", "Divisor could not be 0.");
            dec_rem(ans.base, divr_base);
            if (divr.modulus && divd.sgn != divr.sgn) ans.base = dec_add(ans.sgn, ans.base, false, divr_base, divr.sgn);
        }
        return ans;
        neunet_type_endif
    }
    callback_dec void operator%=(dec_t &&divr) {
        neunet_type_if (neunet_dec_num(dec_t)) *this = *this % net_decimal {divr};
        neunet_type_else *this = *this % divr;
        neunet_type_endif
    }
    callback_dec_arg_ref friend void operator%=(arg &divd, net_decimal &divr) {
        static_assert(std::is_integral_v<arg>, "Dividend should be an integer.");
        int64_t divr_tmp {0};
        if (divr.dec) divr_tmp = dec_to_integer(divr.sgn, divr.base);
        else divr_tmp = dec_to_integer(divr.sgn, divr.frac);
        auto divd_sgn = divd < 0;
        divd %= divr_tmp;
        if (divr.modulus && divd_sgn != divr.sgn) divd += divr_tmp;
    }

    friend std::ostream &operator<<(std::ostream &os, const net_decimal &src) {
        if (src.sgn) os << '-';
        if (src.dec) {
            os << src.base;
            return os;
        }
        auto tmp = dec_div(src.frac.num, src.frac.den, division_precision);
        os << tmp;
        return os;
    }
};

NEUNET_END

_STD_BEGIN



_STD_END
