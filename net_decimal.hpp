NEUNET_BEGIN

struct net_decimal_data final {
    // 456121.484844141
    // {1, 2, 1, 6, 5, 4}, {4, 8, 4, 8, 4, 4, 1, 4, 1}
    net_set<int8_t> it, ft;

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

    net_decimal_data &operator=(const net_decimal_data &src) {
        value_copy(src);
        return *this;
    }

    net_decimal_data &operator=(net_decimal_data &&src) {
        value_move(std::move(src));
        return *this;
    }

    void reset() {
        it.reset();
        ft.reset();
    }

    ~net_decimal_data() { reset(); }
};

bool dec_is_zero(const net_decimal_data &src) { return !(src.it.length || src.ft.length); }

bool dec_is_one(const net_decimal_data &src) { return !src.ft.length && src.it.length == 1 && src.it[0] == 1; }

std::string dec_to_string(bool sgn, const net_decimal_data &src) {
    std::string ans = "";
    if (src.it.length) for (auto i = 0ull; i < src.it.length; ++i) ans.push_back(src.it[src.it.length - 1 - i] + '0');
    else ans.push_back('0');
    if (src.ft.length) {
        ans.push_back('.');
        for (auto i = 0ull; i < src.ft.length; ++i) ans.push_back(src.ft[i] + '0');
    }
    if (sgn && ans != "0") ans = '-' + ans;
    return ans;
}

/* valid digit count = 16
 * value for initialization should be no sign
 */
net_decimal_data dec_init(bool &sgn, long double init_val) {
    uint64_t    it_tmp = std::abs(init_val),
                cnt    = 0;
    long double ft_tmp = std::abs(init_val) - it_tmp;
    auto        p_tmp  = ptr_init<int8_t>(NEUNET_DEC_VLD_DIG);
    net_decimal_data ans;
    sgn = init_val < 0;
    // integer
    while (it_tmp) {
        p_tmp[cnt++] = it_tmp % 10;
        it_tmp      /= 10;
    }
    ans.it = {ptr_sub(cnt, p_tmp, NEUNET_DEC_VLD_DIG, 0ull, cnt - 1), cnt};
    if (cnt >= NEUNET_DEC_VLD_DIG || !ft_tmp) {
        ptr_reset(p_tmp);
        return ans;
    }
    it_tmp = NEUNET_DEC_VLD_DIG - cnt;
    cnt    = 0;
    while (it_tmp--) {
        ft_tmp    *= 10;
        p_tmp[cnt] = ft_tmp;
        ft_tmp    -= p_tmp[cnt++];
    }
    while (!p_tmp[cnt - 1]) --cnt;
    ans.ft = {ptr_sub(cnt, p_tmp, NEUNET_DEC_VLD_DIG, 0ull, cnt - 1), cnt};
    ptr_reset(p_tmp);
    return ans;
}
net_decimal_data dec_init(bool &sgn, const char *init_val) {
    // verify
    auto str_len  = std::strlen(init_val),
         dot_idx  = str_len,
         bgn_idx  = 0ull,
         end_idx  = dot_idx - 1;
    // symbol
    auto dot_flag = false;
    net_decimal_data ans;
    for (auto i = 0ull; i < str_len; ++i) if (init_val[i] < '0' || init_val[i] > '9') {
        if (init_val[i] == '.' && !dot_flag) {
            dot_flag = true;
            dot_idx  = i;
            continue;
        } 
        if (init_val[i] == '-' && !i) {
            sgn = true;
            ++bgn_idx;
            continue;
        }
        if (init_val[i] == '+' && !i) {
            ++bgn_idx;
            continue;
        }
        return ans;
    }
    // zero
    while (end_idx > dot_idx && init_val[end_idx] == '0') --end_idx;
    while (bgn_idx < dot_idx && init_val[bgn_idx] == '0') ++bgn_idx;
    if (end_idx == bgn_idx && dot_idx == bgn_idx) return ans;
    // integer
    auto dig_cnt = dot_idx - bgn_idx;
    if (dig_cnt) {
        ans.it.init(dig_cnt);
        for (auto i = dot_idx; i > bgn_idx; --i) ans.it[dot_idx - i] = init_val[i - 1] - '0';
    }
    if (end_idx < dot_idx) return ans;
    // float
    dig_cnt = end_idx - dot_idx;
    if (dig_cnt) {
        ans.ft.init(dig_cnt);
        for (auto i = dot_idx; i < end_idx; ++i) ans.ft[i - dot_idx] = init_val[i + 1] - '0';
    }
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

// unsigned number digit segment
int8_t dec_add(int8_t fst, int8_t snd, bool &carry) {
    auto ans = fst + snd + carry;
    carry    = false;
    if (ans < 10) return ans;
    carry = true;
    ans  %= 10;
    return ans;
}

int8_t dec_sub(int8_t minu, int8_t subt, bool &carry) {
    subt += carry;
    carry = false;
    if (minu >= subt) return minu - subt;
    carry = true;
    return 10 - subt + minu;
}

int8_t dec_add_sub(int8_t fst, int8_t snd, bool &carry, bool sub = false) { return sub ? dec_sub(fst, snd, carry) : dec_add(fst, snd, carry); }
/* unsigned segment number
 * minuhend (first) segment should be greater than subtrahend (second) segment for true value of parameter subtract
 */
net_decimal_data dec_add_sub(const net_decimal_data &fst, const net_decimal_data &snd, bool sub = false) {
    auto dig_len = std::max(fst.ft.length, snd.ft.length),
         dig_idx = dig_len;
    auto dig_tmp = ptr_init<int8_t>(dig_len);
    auto carry   = false;
    net_decimal_data ans;
    // float
    for (auto i = dig_len; i; --i) {
        auto idx = i - 1;
        auto tmp = dec_add_sub(idx < fst.ft.length ? fst.ft[idx] : 0,
                               idx < snd.ft.length ? snd.ft[idx] : 0,
                               carry, sub);
        if (dig_idx < dig_len || tmp) dig_tmp[--dig_idx] = tmp;
    }
    ans.ft  = {ptr_sub(dig_idx, dig_tmp, dig_len, dig_idx, dig_len - 1), dig_idx};
    // integer
    dig_idx = std::max(fst.it.length, snd.it.length) + 1;
    if (dig_idx > dig_len) {
        ptr_alter(dig_tmp, dig_len, dig_idx, false);
        dig_len = dig_idx;
    } else dig_len = dig_idx;
    dig_idx = 0;
    for (auto i = 0ull; i < dig_len; ++i) dig_tmp[dig_idx++] = dec_add_sub(i < fst.it.length ? fst.it[i] : 0, i < snd.it.length ? snd.it[i] : 0, carry, sub);
    if (carry) dig_tmp[dig_idx++] = carry;
    else while (!dig_tmp[dig_idx - 1]) --dig_idx;
    ans.it  = {ptr_sub(dig_idx, dig_tmp, dig_len, 0ull, dig_idx - 1), dig_idx};
    ptr_reset(dig_tmp);
    return ans;
}

callback_dec_arg int64_t dec_coe_trans(const arg &src, uint64_t ans_len) {
    if constexpr (number_cplx_vld) return src.real() / ans_len + 0.5;
    else return src;
}

uint8_t dec_coe_carry(int64_t &carry) {
    auto tmp = carry;
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

// from low digit
callback_dec_arg net_decimal_data dec_coe_res(const net_set<arg> &src, uint64_t ft_len) {
    auto ans_len = src.length * NEUNET_DEC_DIG_MAX,
         ans_idx = ans_len,
         src_idx = 0ull;
    auto ans_tmp = ptr_init<int8_t>(ans_len);
    auto carry   = 0ll;
    net_decimal_data ans;
    // float
    for (auto i = 0ull; i < ft_len; ++i) {
        if (src_idx < src.length) carry += dec_coe_trans(src[src_idx++], src.length);
        auto tmp = dec_coe_carry(carry);
        if (tmp || ans_idx < ans_len) ans_tmp[--ans_idx] = tmp;
    }
    if (ans_idx < ans_len) ans.ft = {ptr_sub(ans_idx, ans_tmp, ans_len, ans_idx, ans_len - 1), ans_idx};
    ans_idx = 0;
    // integer
    for (auto i = src_idx; i < src.length; ++i) {
        carry += dec_coe_trans(src[i], src.length);
        ans_tmp[ans_idx++] = dec_coe_carry(carry);
    }
    if (carry) while (carry) ans_tmp[ans_idx++] = dec_coe_carry(carry);
    while (!ans_tmp[ans_idx - 1]) --ans_idx;
    if (ans_idx > 0) ans.it = {ptr_sub(ans_idx, ans_tmp, ans_len, 0ull, ans_idx - 1), ans_idx};
    ptr_reset(ans_tmp);
    return ans;
}

// unsigned digit
net_decimal_data dec_mul(const net_decimal_data &fst, const net_decimal_data &snd) {
    auto ft_len  = fst.ft.length + snd.ft.length,
         ans_len = ft_len + fst.it.length + snd.it.length - 1;
    net_set<int64_t> ans_coe(ans_len);
    for (auto i = 0ull; i < fst.ft.length; ++i) {
        auto fst_coe = fst.ft[fst.ft.length - 1 - i];
        for (auto j = 0ull; j < snd.ft.length; ++j) ans_coe[i + j] += fst_coe * snd.ft[snd.ft.length - 1 - j];
        for (auto j = 0ull; j < snd.it.length; ++j) ans_coe[i + j + snd.ft.length] += fst_coe * snd.it[j];
    }
    for (auto i = 0ull; i < fst.it.length; ++i) {
        auto fst_coe = fst.it[i];
        for (auto j = 0ull; j < snd.ft.length; ++j) ans_coe[i + j + fst.ft.length] += fst_coe * snd.ft[snd.ft.length - 1 - j];
        for (auto j = 0ull; j < snd.it.length; ++j) ans_coe[i + j + fst.ft.length + snd.ft.length] += fst_coe * snd.it[j];
    }
    // for (auto tmp : ans_coe) std::cout << tmp << std::endl;
    return dec_coe_res(ans_coe, ft_len);
}

net_set<int8_t> dec_div_coe(const net_decimal_data &src) {
    auto ans = src.it;
    ans.reverse();
    return ans.unit(src.ft);
}

int8_t dec_div(net_set<int8_t> &divd, const net_set<int8_t> &divr) {
    net_set<int8_t> divd_tmp;
    if (divd.length > divr.length) divd_tmp.init(divd.length - 1);
    else divd_tmp.init(divr.length);
    int8_t ans = divd[0] / divr[0];
    // if (!(a[0] % b[0])) --ans;
    for (auto i = 0ull; i < divd.length; ++i) {
        // eliminate the highest order coefficient
        auto tmp_i = i ? i - 1 : i;
        if (i < divr.length) {
            divd_tmp[tmp_i] += divd[i] - ans * divr[i];
            if (i == tmp_i) divd_tmp[tmp_i] *= 10;
        } else divd_tmp[tmp_i] += divd[i];
    }
    for (auto i = divd.length; i < divr.length; ++i) divd_tmp[i - 1] -= ans * divr[i];
    // eliminate
    if (std::abs(divd_tmp[0] / 10) >= std::abs(divr[0])) {
        divd.init(divd_tmp.length + 1, false);
        auto carry   = 0i8;
        auto dvd_dig = divd.length;
        // carrying
        for (auto i = divd_tmp.length; i; --i) {
            --dvd_dig;
            divd[dvd_dig] += carry;
            divd[dvd_dig] += divd_tmp[i - 1] % 10;
            carry          = divd_tmp[i - 1] / 10;
        }
        divd[0] = carry;
        ans    += dec_div(divd, divr);
    } else divd = std::move(divd_tmp); 
    return ans;
}
net_decimal_data dec_div(const net_decimal_data &divd, const net_decimal_data &divr, uint64_t prec) {
    auto divd_coe = dec_div_coe(divd),
         divr_coe = dec_div_coe(divr);
    auto it_seg_c = 1ull;
    if (divd.it.length >= divr.it.length) it_seg_c = divd.it.length - divr.it.length + 1;
    auto idx_tool = 0ull,
         ft_len   = prec + 2;
    net_set<int8_t> ans_set(it_seg_c + ft_len);
    if (divd.it.length < divr.it.length) idx_tool = divr.it.length - divd.it.length;
    idx_tool = ans_set.length - idx_tool;
    for (auto i = idx_tool; i > 2; --i) ans_set[i - 1] = dec_div(divd_coe, divr_coe);
    return dec_coe_res(ans_set, ft_len);
}

bool dec_mod_verify(net_set<int8_t> &divd, const net_set<int8_t> &divr) {
    int64_t carry = 0ull;
    for (auto i = divd.length; i > 1; --i) {
        auto idx  = i - 1;
        carry    += divd[idx];
        divd[idx] = dec_coe_carry(carry);
    }
    divd[0] += carry;
    if (divd[0] < 0) return false;
    for (auto i = 0ull; i < divd.length; ++i) {
        if (divd[i] > divr[i]) return false;
        if (divd[i] < divr[i]) return true;
    }
    return true;
}

int8_t dec_mod(net_set<int8_t> &divd, const net_set<int8_t> &divr, int8_t coe = 0) {
    if (divd.length != divr.length || dec_mod_verify(divd, divr)) return 0;
    auto tmp = divd[0];
    auto ans = tmp / divr[0];
    if ((-tmp) == coe) ans /= 2;
    if (!ans) --ans;
    for (auto i = 0ull; i < divd.length; ++i) divd[i] -= ans * divr[i];
    return ans + dec_mod(divd, divr, tmp);
}
net_decimal_data dec_mod(net_decimal_data &divd_rem, const net_decimal_data &divr) {
    if (divd_rem.ft.length || divr.ft.length) return {};
    auto comp_res = dec_comp(divd_rem, divr);
    auto sgn_tmp  = false;
    if (comp_res == NEUNET_DEC_CMP_LES) return dec_init(sgn_tmp, 0.);
    if (comp_res == NEUNET_DEC_CMP_EQL) return dec_init(sgn_tmp, 1);
    auto divd_coe = dec_div_coe(divd_rem);
    auto divr_coe = dec_div_coe(divr);
    auto dig_cnt  = divd_rem.it.length - divr.it.length + 1,
         divd_idx = 0ull;
    net_set<int64_t> ans_set(dig_cnt); // --dig_cnt;
    while (divd_coe.length > divr_coe.length) ans_set[--dig_cnt] = dec_div(divd_coe, divr_coe);
    // last digit
    ans_set[--dig_cnt] = dec_mod(divd_coe, divr_coe);
    // for (auto tmp : ans_set) std::cout << (int)tmp << std::endl;
    // std::cout << "------" << std::endl;
    // for (auto tmp : divd_coe) std::cout << (int)tmp << std::endl;
    divd_coe.reverse();
    divd_rem = dec_coe_res(divd_coe, divd_rem.ft.length);
    return dec_coe_res(ans_set, 0);
}

net_set<uint64_t> dec_fft_rev(uint64_t fst_len, uint64_t snd_len) {
    net_set<uint64_t> ans(1ull << num_bit_cnt((fst_len + snd_len) - 1));
    for (int i = 0; i < ans.length; ++i) ans[i] = (ans[i >> 1] >> 1) | ((i & 1) << (std::max(uint64_t(std::log2(ans.length - 1) + 0.5), 1ull) - 1));
    return ans;
}

void dec_fft_coe(net_set<std::complex<long double>> &fst_coe, net_set<std::complex<long double>> &snd_coe, const net_decimal_data &fst, const net_decimal_data &snd, const net_set<uint64_t> &rev, uint64_t fst_len, uint64_t snd_len) {
    fst_coe.init(rev.length, false);
    snd_coe.init(rev.length, false);
    for (auto i = 0ull; i < rev.length; ++i) {
        auto idx = rev[i];
        if (idx < fst_len) {
            if (idx < fst.ft.length) fst_coe[i] = fst.ft[fst.ft.length - idx - 1];
            else fst_coe[i] = fst.it[idx - fst.ft.length];
        }
        if (idx < snd_len) {
            if (idx < snd.ft.length) snd_coe[i] = snd.ft[snd.ft.length - idx - 1];
            else snd_coe[i] = snd.it[idx - snd.ft.length];
        }
    }
}
void dec_fft_coe(net_set<std::complex<long double>> &coe, const net_set<uint64_t> &rev) { for (auto i = 0ull; i < rev.length; ++i) if (i < rev[i]) std::swap(coe[i], coe[rev[i]]); }

void dec_fft(net_set<std::complex<long double>> &rev_coe, bool inv = false) {
    auto half = inv ? -1 : 1;
    for (auto i = 1ull; i < rev_coe.length; i <<= 1) {
        std::complex<long double> w{std::cos(NEUNET_DEC_PI / i),
                                    std::sin(NEUNET_DEC_PI / i) * half};
        for (auto j = 0ull; j < rev_coe.length; j += (i << 1)) {
            std::complex<long double> w_x{1, 0};
            for (auto k = 0ull; k < i; ++k, w_x *= w) {
                auto x = rev_coe[j + k],
                     y = w_x * rev_coe[i + j + k];
                rev_coe[j + k]     = x + y;
                rev_coe[i + j + k] = x - y;
            }
        }
    }
}

net_decimal_data dec_fft_mul(const net_decimal_data &fst, const net_decimal_data &snd) {
    auto ft_len  = fst.ft.length + snd.ft.length,
         fst_len = fst.ft.length + fst.it.length,
         snd_len = snd.ft.length + snd.it.length;
    auto coe_rev = dec_fft_rev(fst_len, snd_len);
    net_set<std::complex<long double>> fst_coe, snd_coe, ans_coe(coe_rev.length);
    dec_fft_coe(fst_coe, snd_coe, fst, snd, coe_rev, fst_len, snd_len);
    dec_fft(fst_coe);
    dec_fft(snd_coe);
    for (auto i = 0ull; i < coe_rev.length; ++i) ans_coe[i] = fst_coe[i] * snd_coe[i];
    dec_fft_coe(ans_coe, coe_rev);
    dec_fft(ans_coe, true);
    // for (auto tmp : ans_coe) std::cout << (int)(tmp.real() / ans_coe.length + 0.5) << std::endl;
    return dec_coe_res(ans_coe, ft_len);
}

void dec_truncate(net_decimal_data &src, uint64_t acc) {
    if (acc >= src.ft.length) return;
    if (acc + 1 == src.ft.length && src.ft[acc] < 5) return src.ft.init(acc);
    acc += 2;
    if (src.ft[acc - 1] >= 5) {
        net_decimal_data carry;
        carry.ft.init(acc);
        carry.ft[acc - 1] = 5;
        src = dec_add_sub(src, carry);
    }
    acc -= 2;
    if (src.ft[acc] < 5) src.ft.init(acc);
    else src.ft.init(acc + 1);
}

bool dec_gcd(net_decimal_data &fst, net_decimal_data &snd) {
    auto comp = dec_comp(fst, snd);
    if (comp == NEUNET_DEC_CMP_EQL) return false;
    if (comp == NEUNET_DEC_CMP_GTR) {
        dec_mod(fst, snd);
        // std::cout << dec_to_string(false, fst) << '\n' << dec_to_string(false, snd) << '\n' << std::endl;
        if (dec_is_zero(fst)) return true;
    }
    if (comp == NEUNET_DEC_CMP_LES) {
        dec_mod(snd, fst);
        // std::cout << dec_to_string(false, fst) << '\n' << dec_to_string(false, snd) << '\n' << std::endl;
        if (dec_is_zero(snd)) return false;
    }
    return dec_gcd(fst, snd);
}

// unsigned value, unsafe
class net_decimal_base final {
public:
    net_decimal_base() {}
    net_decimal_base(const net_decimal_base &src) { val = src.val; }
    net_decimal_base(net_decimal_base &&src) { val = std::move(src.val); }

    bool init(const std::string &src) {
        auto sgn = false;
        val = dec_init(sgn, src.c_str());
        return sgn;
    }
    bool init(long double src) {
        auto sgn = false;
        val = dec_init(sgn, src);
        return sgn;
    }

    bool is_zero() const { return dec_is_zero(val); }

    bool is_one() const { return dec_is_one(val); }

    net_decimal_base integer_part() const {
        net_decimal_base ans;
        if (val.it.length) ans.val.it = val.it;
        return ans;
    }

    net_decimal_base float_part() const {
        net_decimal_base ans;
        if (val.ft.length) ans.val.ft = val.ft;
        return ans;
    }

    std::string to_string(bool sgn) const { return dec_to_string(sgn, val); }

    uint64_t to_integer() const {
        if (!val.it.length) return 0;
        auto ans = 0ull;
        for (auto i = val.it.length; i; --i) {
            ans *= 10;
            ans += val.it[i - 1];
        } 
        return ans;
    }

    long double to_float() const {
        long double ans_it = to_integer(),
                    ans_ft = 0;
        for (auto i = 0ull; i < val.ft.length; ++i) {
            ans_ft += val.ft[i];
            ans_ft /= 10;
        }
        return ans_it + ans_ft;
    }

    int comp(const net_decimal_base &src) const { return dec_comp(val, src.val); }

    net_decimal_base gcd(const net_decimal_base &src) const {
        net_decimal_base ans;
        if (is_one() || src.is_one()) {
            ans.init(1);
            return ans;
        }
        auto fst_in = val,
             snd_in = src.val;
        auto ans_id = dec_gcd(fst_in, snd_in);
        if (ans_id) ans.val = std::move(snd_in);
        else ans.val = std::move(fst_in);
        return ans;
    }

    net_decimal_base mod(net_decimal_base &div_rem, const net_decimal_base &divr) const {
        net_decimal_base ans;
        div_rem.val = val;
        ans.val     = dec_mod(div_rem.val, divr.val);
        return ans;
    }

    net_decimal_base exp10(net_decimal_base &e10) const {
        e10.val.ft.reset();
        e10.val.it.init(val.ft.length + 1);
        e10.val.it[val.ft.length] = 1;
        net_decimal_base ans;
        ans.val.it = val.ft;
        ans.val.it.reverse();
        ans.val.it = ans.val.it.unit(val.it);
        return ans;
    }

    uint64_t float_digit_cnt() const { return val.ft.length; }

    uint64_t integer_digit_cnt() const { return val.it.length; }

    void reset() { val.reset(); }

    ~net_decimal_base() { reset(); }

    friend net_decimal_base operator+(const net_decimal_base &fst, const net_decimal_base &snd) {
        net_decimal_base ans;
        ans.val = dec_add_sub(fst.val, snd.val);
        return ans;
    }

    friend net_decimal_base operator-(const net_decimal_base &fst, const net_decimal_base &snd) {
        net_decimal_base ans;
        ans.val = dec_add_sub(fst.val, snd.val, true);
        return ans;
    }

    friend net_decimal_base operator*(const net_decimal_base &fst, const net_decimal_base &snd) {
        net_decimal_base ans;
        if (base_fft) ans.val = dec_fft_mul(fst.val, snd.val);
        else ans.val = dec_mul(fst.val, snd.val);
        return ans;
    }

    friend net_decimal_base operator/(const net_decimal_base &fst, const net_decimal_base &snd) {
        net_decimal_base ans;
        ans.val = dec_div(fst.val, snd.val, base_acc);
        return ans;
    }

    net_decimal_base &operator=(const net_decimal_base &src) {
        val = src.val;
        return *this;
    }
    net_decimal_base &operator=(net_decimal_base &&src) {
        val = std::move(src.val);
        return *this;
    }

private: net_decimal_data val;

public:
    static std::atomic_bool base_fft;

    static std::atomic_uint64_t base_acc;
};

// unsigned value, unsafe
// The denominator and fractor should be integer in fraction format
class net_decimal_frac final {
private:
    void value_copy(const net_decimal_frac &src) {
        num = src.num;
        den = src.den;
    }

    void value_move(net_decimal_frac &&src) {
        num = std::move(src.num);
        den = std::move(src.den);
    }

public:
    net_decimal_frac() {}
    net_decimal_frac(const net_decimal_frac &src) { value_copy(src); }
    net_decimal_frac(net_decimal_frac &&src) { value_move(std::move(src)); }

    void init(const net_decimal_base &src) {
        num = src.exp10(den);
        reduct();
    }
    bool init(const net_decimal_base &src_num, const net_decimal_base &src_den) {
        if (src_den.is_zero() || src_num.float_digit_cnt() || src_den.float_digit_cnt()) return false;
        num = src_num;
        den = src_den;
        return true;
    }

    bool is_valid() const { return !den.is_zero(); }

    bool is_zero() {
        if (!num.is_zero()) return false;
        den.init(1);
        return true;
    }

    bool is_one() {
        if (num.comp(den) != NEUNET_DEC_CMP_EQL) return false;
        num.init(1);
        den.init(1);
        return true;
    }

    void reduct() {
        auto gcd_val = den.gcd(num);
        num = num / gcd_val;
        den = den / gcd_val;
    }

    int comp(const net_decimal_frac &src) {
        auto fst = num * src.den,
             snd = src.num * den;
        return fst.comp(snd);
    }

    std::string to_string(bool sgn) const { return (num / den).to_string(sgn); }

    uint64_t to_integer() const { return (num / den).to_integer(); }

    long double to_float() const { return (num / den).to_float(); }

    // return integer part and the output parameter of numerator
    net_decimal_frac integer_part(net_decimal_frac &float_part) const {
        net_decimal_frac ans;
        ans.den.init(1);
        float_part.den = den;
        ans.num = num.mod(float_part.num, den);
        return ans;
    }

    bool is_integer() {
        reduct();
        return den.is_one();
    }

    const net_decimal_base &numerator() const { return num; }

    const net_decimal_base &denominator() const { return den; }

    void reset() {
        num.reset();
        den.reset();
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

    friend net_decimal_frac operator+(net_decimal_frac &fst, net_decimal_frac &snd) {
        net_decimal_frac ans;
        if (!(fst.is_valid() && snd.is_valid())) return ans;
        ans.den = fst.den * snd.den;
        ans.num = fst.num * snd.den + snd.num * fst.den;
        return ans;
    }

    friend net_decimal_frac operator-(net_decimal_frac &fst, net_decimal_frac &snd) {
        net_decimal_frac ans;
        if (!(fst.is_valid() && snd.is_valid())) return ans;
        ans.den = fst.den * snd.den;
        ans.num = fst.num * snd.den - snd.num * fst.den;
        return ans;
    }

    friend net_decimal_frac operator*(net_decimal_frac &fst, net_decimal_frac &snd) {
        net_decimal_frac ans;
        if (!(fst.is_valid() && snd.is_valid())) return ans;
        ans.num = fst.num * snd.num;
        ans.den = fst.den * snd.den;
        return ans;
    }

    friend net_decimal_frac operator/(net_decimal_frac &fst, net_decimal_frac &snd) {
        net_decimal_frac ans;
        if (!(fst.is_valid() && snd.is_valid())) return ans;
        ans.num = fst.num * snd.den;
        ans.den = fst.den * snd.num;
        return ans;
    }

private: net_decimal_base num, den; // numerator / denominator
};

callback_dec_ins ins dec_add(bool &ans_sgn, ins &fst, bool fst_sgn, ins &snd, bool snd_sgn) {
    if (fst.is_zero()) {
        ans_sgn = snd_sgn;
        return snd;
    }
    if (snd.is_zero()) {
        ans_sgn = fst_sgn;
        return fst;
    }
    if (fst_sgn != snd_sgn) {
        auto comp_ans = snd.comp(fst);
        if (comp_ans == NEUNET_DEC_CMP_GTR) {
            if (snd_sgn) ans_sgn = true;
            return snd - fst;
        }
        if (comp_ans == NEUNET_DEC_CMP_LES) {
            if (fst_sgn) ans_sgn = true;
            return fst - snd;
        }
    }
    if (fst_sgn) ans_sgn = true;
    return fst + snd;
}

callback_dec_ins ins dec_sub(bool &ans_sgn, ins &minu, bool minu_sgn, ins &subt, bool subt_sgn) {
    if (subt.is_zero()) {
        ans_sgn = minu_sgn;
        return minu;
    }
    if (minu.is_zero()) {
        ans_sgn = !subt_sgn;
        return subt;
    }
    if (minu_sgn == subt_sgn) {
        auto comp_ans = minu.comp(subt);
        if (comp_ans == NEUNET_DEC_CMP_GTR) {
            if (minu_sgn) ans_sgn = true;
            return minu - subt;
        }
        if (comp_ans == NEUNET_DEC_CMP_LES) {
            if (!minu_sgn) ans_sgn = true;
            return subt - minu;
        }
    }
    if (minu_sgn) ans_sgn = true;
    return minu + subt;
}

callback_dec_ins ins dec_mul(bool &ans_sgn, ins &fst, bool fst_sgn, ins &snd, bool snd_sgn) {
    ans_sgn = false;
    if (fst.is_zero() || snd.is_zero()) return {};
    ans_sgn = snd_sgn != fst_sgn;
    if (fst.is_one()) return snd;
    if (snd.is_one()) return fst;
    return fst * snd;
}

class net_decimal {
public:
    bool modulus = false;

    static std::atomic_bool fft_mode;

    static std::atomic_uint64_t accuracy;

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

    bool is_zero() { return dec ? base.is_zero() : frac.is_zero(); }

    bool is_one() { return dec ? base.is_one() : frac.is_one(); }

    // ln(x + 1) = (-1 < x <= 1)
    net_decimal ln_2() {
        net_decimal b {1},
                    p {0},
                    q {b},
                    x {*this - b},
                    o {x},
                    c {b},
                    r {p};
        std::string m {"0"},
                    n {"1"};
        while (m != n) {
            // std::cout << n << std::endl;
            m = std::move(n);
            if (b == q) p += x;
            else {
                p  = b * p + c * x * q;
                q *= b;
            }
            c.sgn = !c.sgn;
            x    *= o;
            r     = p / q;
            n     = r.to_string();
            ++b;
        }
        return r;
    }
    static net_decimal ln_4() {
        net_decimal b {1},
                    c {2},
                    s {"0.6"},
                    o {"0.36"},
                    p {0},
                    q {b},
                    r {p};
        std::string m {"0"},
                    n {"1"};
        while (m != n) {
            // std::cout << n << std::endl;
            m = std::move(n);
            if (b == q) p += s;
            else {
                p  = b * p + s * q;
                q *= b;
            }
            s *= o;
            b += 2;
            r  = p / q;
            n  = r.to_string();
        }
        return 2 * r;
    }
    net_decimal dec_ln_proto() {
        net_decimal b {1},
                    p {0},
                    q {b},
                    o {p},
                    u {*this - 1},
                    v {*this + 1},
                    h {u * u},
                    s {v * v},
                    r {p};
        std::string m {"0"},
                    n {"1"};
        while (m != n) {
            // std::cout << n << std::endl;
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
            r  = p / q;
            n  = r.to_string();
        }
        return r * 2;
    }

public:
    net_decimal() {}
    callback_arg net_decimal(const arg &src) {
        static_assert(std::is_arithmetic_v<arg> ||
                      std::is_same_v<arg, arr_t_v(char, arg)> ||
                      std::is_same_v<std::string, arg>,
                      "A arithmetic or string type value is needed.");
        sgn = base.init(src);
    }
    net_decimal(const net_decimal &src) { value_copy(src); }
    net_decimal(net_decimal &&src) { value_move(std::move(src)); }

    std::string to_string() const { return dec ? base.to_string(sgn) : frac.to_string(sgn); }

    uint64_t to_integer() const { return dec ? base.to_integer() : frac.to_integer(); }

    long double to_float() const {
        auto ans = dec ? base.to_float() : frac.to_float();
        if (sgn) return -ans;
        return ans;
    }

    net_decimal integer_part() const {
        net_decimal ans;
        if (dec) ans.base = base.integer_part();
        else {
            ans.dec  = false;
            ans.frac = frac.integer_part(ans.frac);
        }
        ans.sgn = sgn;
        return ans;
    }

    net_decimal float_part() const {
        net_decimal ans;
        if (dec) ans.base = base.float_part();
        else {
            ans.dec = false;
            frac.integer_part(ans.frac);
        }
        ans.sgn = sgn;
        return ans;
    }

    int compare(net_decimal &src) {
        if (src.dec == dec && dec) return base.comp(src.base);
        if (dec) if (!frac.is_valid()) frac.init(base);
        if (src.dec) if (!src.frac.is_valid()) src.frac.init(src.base);
        return frac.comp(src.frac);
    }

    net_decimal ln() {
        if (is_zero() || is_one() || sgn) return 0.;
        net_decimal ln4  = 0.,
                    base = *this;
        while (base > 1) {
            base *= 0.25;
            ++ln4;
        }
        ln4 *= ln_4();
        return base.ln_2() + ln4;
    }

    void reset() {
        sgn = false;
        dec = true;
        base.reset();
        frac.reset();
    }

    ~net_decimal() { reset(); }

protected:
    bool sgn = false,
         dec = true;

    net_decimal_base base;

    net_decimal_frac frac;

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
        // first is arithmetic type
        dec_num_if (fst_dec_t)
        return net_decimal {fst} <=> snd;
        // second is arithmetic type
        dec_num_elif (snd_dec_t)
        return fst <=> net_decimal {snd};
        // both are decimal
        dec_num_else
        if (fst.sgn != snd.sgn) return fst.sgn ? std::strong_ordering::less : std::strong_ordering::greater;
        auto cmp_res = fst.compare(snd);
        if (cmp_res == NEUNET_DEC_CMP_LES) return fst.sgn ? std::strong_ordering::greater : std::strong_ordering::less;
        if (cmp_res == NEUNET_DEC_CMP_GTR) return fst.sgn ? std::strong_ordering::less : std::strong_ordering::greater;
        return std::strong_ordering::equal;
        dec_num_endif
    }
    callback_dec_s friend bool operator==(fst_dec_t &&fst, snd_dec_t &&snd) {
        // first is arithmetic type
        dec_num_if (fst_dec_t)
        return net_decimal {fst} == snd;
        // second is arithmetic type
        dec_num_elif (snd_dec_t)
        return fst == net_decimal {snd};
        // both are decimal
        dec_num_else
        return fst.compare(snd) == NEUNET_DEC_CMP_EQL;
        dec_num_endif
    }

    /* TODO fraction has the top privilege in accuracy limit */

    net_decimal operator+() { return *this; }
    callback_dec_s friend net_decimal operator+(fst_dec_t &&fst, snd_dec_t &&snd) {
        // first is arithmetic type
        dec_num_if (fst_dec_t)
        return net_decimal {fst} + snd;
        // second is arithmetic type
        dec_num_elif (snd_dec_t)
        return fst + net_decimal {snd};
        // both are decimal
        dec_num_else
        net_decimal ans;
        if (fst.dec && snd.dec) {
            ans.base = dec_add(ans.sgn, fst.base, fst.sgn, snd.base, snd.sgn);
            return ans;
        }
        if (fst.dec) fst.frac.init(fst.base);
        if (snd.dec) snd.frac.init(snd.base);
        ans.frac = dec_add(ans.sgn, fst.frac, fst.sgn, snd.frac, snd.sgn);
        ans.dec  = false;
        return ans;
        dec_num_endif
    }
    callback_dec void operator+=(dec_t &&src) {
        dec_num_if (dec_t) *this = *this + net_decimal {src};
        dec_num_else *this = *this + src;
        dec_num_endif
    }
    callback_dec_arg friend void operator+=(arg &fst, const net_decimal &snd) { fst += snd.to_float(); }
    net_decimal &operator++() {
        *this += net_decimal {1};
        return *this;
    }
    net_decimal operator++(int) {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    net_decimal operator-() { return net_decimal {} - *this; }
    callback_dec_s friend net_decimal operator-(fst_dec_t &&minu, snd_dec_t &&subt) {
        // first is arithmetic type
        dec_num_if (fst_dec_t)
        return net_decimal {minu} - subt;
        // second is arithmetic type
        dec_num_elif (snd_dec_t)
        return minu - net_decimal {subt};
        // both are decimal
        dec_num_else
        net_decimal ans;
        if (minu.dec && subt.dec) {
            ans.base = dec_sub(ans.sgn, minu.base, minu.sgn, subt.base, subt.sgn);
            return ans;
        }
        if (minu.dec) minu.frac.init(minu.base);
        if (subt.dec) subt.frac.init(subt.base);
        ans.frac = dec_sub(ans.sgn, minu.frac, minu.sgn, subt.frac, subt.sgn);
        ans.dec  = false;
        return ans;
        dec_num_endif
    }
    callback_dec void operator-=(dec_t &&src) {
        dec_num_if (dec_t) *this = *this - net_decimal {src};
        dec_num_else *this = *this - src;
        dec_num_endif
    }
    callback_dec_arg friend void operator-=(arg &fst, const net_decimal &snd) { fst -= snd.to_float(); }
    net_decimal &operator--() {
        *this -= net_decimal {1};
        return *this;
    }
    net_decimal operator--(int) {
        auto tmp = *this;
        --(*this);
        return tmp;
    }

    callback_dec_s friend net_decimal operator*(fst_dec_t &&fst, snd_dec_t &&snd) {
        // first is arithmetic type
        dec_num_if (fst_dec_t)
        return net_decimal {fst} * snd;
        // second is arithmetic type
        dec_num_elif (snd_dec_t)
        return fst * net_decimal {snd};
        // both are decimal
        dec_num_else
        net_decimal ans;
        if (fst.dec && snd.dec) {
            ans.base = dec_mul(ans.sgn, fst.base, fst.sgn, snd.base, snd.sgn);
            return ans;
        }
        if (fst.dec) fst.frac.init(fst.base);
        if (snd.dec) snd.frac.init(snd.base);
        ans.frac = dec_mul(ans.sgn, fst.frac, fst.sgn, snd.frac, snd.sgn);
        ans.dec  = false;
        return ans;
        dec_num_endif
    }
    callback_dec void operator*=(dec_t &&src) {
        dec_num_if (dec_t) *this = *this * net_decimal {src};
        dec_num_else *this = *this * src;
        dec_num_endif
    }
    callback_dec_arg friend void operator*=(arg &fst, const net_decimal &snd) { fst *= snd.to_float(); }

    callback_dec_s friend net_decimal operator/(fst_dec_t &&divd, snd_dec_t &&divr) {
        // first is arithmetic type
        dec_num_if (fst_dec_t)
        return net_decimal {divd} / divr;
        // second is arithmetic type
        dec_num_elif (snd_dec_t)
        return divd / net_decimal {divr};
        // both are decimal
        dec_num_else
        if (divr.is_zero()) {
            std::cerr << "Divisor could not be 0." << std::endl;
            std::abort();
        }
        if (divd.is_zero()) return {};
        net_decimal ans;
        auto ans_sgn = divr.sgn != divd.sgn;
        if (divr.is_one()) {
            ans     = divd;
            ans.sgn = ans_sgn;
            return ans;
        }
        if (!divd.frac.is_valid()) divd.frac.init(divd.base);
        if (!divr.frac.is_valid()) divr.frac.init(divr.base);
        ans.dec  = false;
        ans.sgn  = ans_sgn;
        ans.frac = divd.frac / divr.frac;
        return ans;
        dec_num_endif
    }
    callback_dec void operator/=(dec_t &&divr) {
        dec_num_if (dec_t) *this = *this / net_decimal {divr};
        dec_num_else *this = *this / divr;
        dec_num_endif
    }
    callback_dec_arg friend void operator/=(arg &fst, const net_decimal &snd) { fst /= snd.to_float(); }

    callback_dec_s friend net_decimal operator%(fst_dec_t &&divd, snd_dec_t &&divr) {
        // first is arithmetic type
        dec_num_if (fst_dec_t)
        return net_decimal {divd} / divr;
        // second is arithmetic type
        dec_num_elif (snd_dec_t)
        return divd / net_decimal {divr};
        // both are decimal
        dec_num_else
        if ((divd.dec ? divd.base.float_digit_cnt() : !divd.frac.is_integer()) ||
            (divr.dec ? divr.base.float_digit_cnt() : !divr.frac.is_integer())) {
            std::cerr << "Dividend and divisor should be integer." << std::endl;
            std::abort();
        }
        net_decimal ans;
        ans.sgn = divd.sgn != divr.sgn;
        if (divd.dec && divr.dec) {
            divd.base.mod(ans.base, divr.base);
            if (divr.modulus && ans.sgn) return ans + divr;
            return ans;
        }
        if (divd.dec) divd.base.mod(ans.base, divr.frac.numerator());
        if (divr.dec) divd.frac.numerator().mod(ans.base, divr.base);
        if (divr.modulus && ans.sgn) return ans + divr;
        return ans;
        dec_num_endif
    }
    callback_dec void operator%=(dec_t &&divr) {
        dec_num_if (dec_t) *this = *this % net_decimal {divr};
        dec_num_else *this = *this % divr;
        dec_num_endif
    }
    callback_dec_arg friend void operator%=(arg &fst, net_decimal &divr) {
        if ((divr.dec ? divr.base.float_digit_cnt() : !divr.frac.is_integer()) || !std::is_integral_v<arg>) {
            std::cerr << "Dividend and divisor should be integer." << std::endl;
            std::abort();
        }
        auto sgn_coe  = 1;
        auto divd_sgn = fst < 0;
        if (divd_sgn != divr.sgn) sgn_coe = -1;
        fst %= divr.to_integer();
        fst *= sgn_coe;
        if (divr.modulus && fst < 0) fst += divr;
    }

    friend std::ostream &operator<<(std::ostream &os, const net_decimal &src) {
        os << src.to_string();
        return os;
    }
};

NEUNET_END

_STD_BEGIN



_STD_END
