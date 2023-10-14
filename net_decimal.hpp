NEUNET_BEGIN

uint64_t dec_add_prec(uint64_t fst_prec, uint64_t snd_prec) {
    if (fst_prec && snd_prec) return fst_prec > snd_prec ? fst_prec : snd_prec;
    if (fst_prec) return fst_prec;
    return snd_prec;
}

// valid digit count = 16
void dec_init(uint64_t &precision, bool &sgn, net_set<uint64_t> &it, net_set<uint64_t> &ft, const long double init_val) {
    sgn = init_val < 0;
    uint64_t    it_tmp = std::abs(init_val),
                cnt    = 0;
    long double ft_tmp = std::abs(init_val) - it_tmp;
    if (it_tmp) {
        it.init(1, false);
        it[0] = it_tmp;
    }
    uint64_t it_seg_cnt = 0,
             ft_dig_cnt = 0;
    while (it_tmp) {
        ++it_seg_cnt;
        it_tmp /= 10;
    }
    precision = 0;
    if (it_seg_cnt >= NEUNET_DEC_VLD_DIG || !ft_tmp) return;
    ft_dig_cnt = NEUNET_DEC_VLD_DIG - it_seg_cnt;
    precision  = ft_dig_cnt;
    ft.init(1, false);
    while (ft_dig_cnt--) {
        ft_tmp *= 10;
        ft[0]  *= 10;
        ft[0]  += ft_tmp;
        ft_tmp -= (int)ft_tmp;
        ++cnt;
    }
    while (cnt++ < NEUNET_DEC_DIG_MAX) ft[0] *= 10;
}

bool dec_init(uint64_t &precision, bool &sgn, net_set<uint64_t> &it, net_set<uint64_t> &ft, const char *init_val) {
    // verify
    auto str_tmp  = std::string(init_val);
    auto dot_idx  = str_tmp.length(),
         bgn_idx  = 0ull,
         end_idx  = dot_idx - 1;
    // symbol
    auto dot_flag = false;
    for (auto i = 0ull; i < str_tmp.length(); ++i) if (str_tmp[i] < '0' || str_tmp[i] > '9') {
        if (str_tmp[i] == '.' && !dot_flag) {
            dot_flag = true;
            dot_idx  = i;
            continue;
        } 
        if (str_tmp[i] == '-' && !i) {
            sgn = true;
            ++bgn_idx;
            continue;
        }
        if (str_tmp[i] == '+' && !i) {
            ++bgn_idx;
            continue;
        }
        return false;
    }
    // zero
    while (end_idx > dot_idx && str_tmp[end_idx] == '0') --end_idx;
    while (bgn_idx < dot_idx && str_tmp[bgn_idx] == '0') ++bgn_idx;
    if (end_idx == bgn_idx && dot_idx == bgn_idx) return true;
    uint64_t seg_cnt = 0,
             seg_tmp = 0,
             tmp_cnt = 0;
    // integer
    if (bgn_idx < dot_idx) {
        tmp_cnt = dot_idx - bgn_idx;
        seg_cnt = tmp_cnt / NEUNET_DEC_DIG_MAX;
        if (tmp_cnt % NEUNET_DEC_DIG_MAX) ++seg_cnt;
        it.init(seg_cnt, false);
        seg_cnt = 0;
        tmp_cnt = 1;
        for (auto i = dot_idx; i > bgn_idx; --i) {
            seg_tmp += (str_tmp[i - 1] - '0') * tmp_cnt;
            tmp_cnt *= 10;
            if (tmp_cnt == NEUNET_DEC_SEG_MAX) {
                it[seg_cnt++] = seg_tmp;
                tmp_cnt       = 1;
                seg_tmp       = 0;
            }
        }
        if (seg_tmp) {
            it[seg_cnt] = seg_tmp;
            seg_tmp     = 0;
        }
        seg_cnt = 0;
    }
    // float
    if (end_idx > dot_idx) {
        tmp_cnt   = end_idx - dot_idx;
        precision = tmp_cnt;
        seg_cnt   = tmp_cnt / NEUNET_DEC_DIG_MAX;
        if (tmp_cnt % NEUNET_DEC_DIG_MAX) ++seg_cnt;
        ft.init(seg_cnt, false);
        seg_cnt = 0;
        tmp_cnt = 0;
        for (auto i = dot_idx + 1; i <= end_idx; ++i) {
            seg_tmp *= 10;
            seg_tmp += str_tmp[i] - '0';
            ++tmp_cnt;
            if (tmp_cnt == NEUNET_DEC_DIG_MAX) {
                ft[seg_cnt++] = seg_tmp;
                tmp_cnt       = 0;
                seg_tmp       = 0;
            }
        }
        if (seg_tmp) {
            while (tmp_cnt++ < NEUNET_DEC_DIG_MAX) seg_tmp *= 10;
            ft[seg_cnt] = seg_tmp;
        }
    }
    return true;
}

/* unsigned decimal digit part
 * NEUNET_DEC_CMP_EQL first is equal to second
 * NEUNET_DEC_CMP_LES first is less than second
 * NEUNET_DEC_CMP_GTR first is greater than second
 */
int dec_compare(const net_set<uint64_t> &fst_it, const net_set<uint64_t> &fst_ft, const net_set<uint64_t> &snd_it, const net_set<uint64_t> &snd_ft) {
    if (fst_it.length > snd_it.length) return NEUNET_DEC_CMP_GTR;
    if (fst_it.length < snd_it.length) return NEUNET_DEC_CMP_LES;
    for (auto i = fst_it.length; i; --i) {
        auto idx = i - 1;
        if (fst_it[idx] > snd_it[idx]) return NEUNET_DEC_CMP_GTR;
        if (fst_it[idx] < snd_it[idx]) return NEUNET_DEC_CMP_LES;
    }
    for (auto i = 0; i < fst_ft.length; ++i) {
        if (i == snd_ft.length) return NEUNET_DEC_CMP_GTR;
        if (fst_ft[i] > snd_ft[i]) return NEUNET_DEC_CMP_GTR;
        if (fst_ft[i] < snd_ft[i]) return NEUNET_DEC_CMP_LES;
    }
    if (fst_ft.length == snd_ft.length) return NEUNET_DEC_CMP_EQL;
    else return NEUNET_DEC_CMP_LES;
}

// unsigned number digit segment
uint64_t dec_add(bool &carry, uint64_t fst, uint64_t snd) {
    auto dif = NEUNET_DEC_SEG_MAX - fst,
         ans = 0ull;
    auto p_c = carry;
    carry    = snd >= dif;
    if (carry) {
        ans = snd - dif;
        if (p_c) ++ans;
    } else {
        ans = fst + snd;
        if (p_c) ans = dec_add(carry, ans, p_c);
    }
    return ans;
}

uint64_t dec_sub(bool &carry, uint64_t minu, uint64_t subt) {
    auto p_c = carry;
    auto ans = 0ull;
    carry    = subt > minu;
    if (carry) {
        ans = NEUNET_DEC_SEG_MAX - subt + minu;
        if (p_c) --ans;
    } else {
        ans = minu - subt;
        if (p_c) ans = dec_sub(carry, ans, p_c);
    }
    return ans;
}

uint64_t dec_add_sub(bool &carry, uint64_t fst, uint64_t snd, bool subtract) { return subtract ? dec_sub(carry, fst, snd) : dec_add(carry, fst, snd); }
/* unsigned segment number
 * minuhend (first) segment should be greater than subtrahend (second) segment for true value of parameter subtract
 */
void dec_add_sub(net_set<uint64_t> &ans_it, net_set<uint64_t> &ans_ft, const net_set<uint64_t> &fst_it, const net_set<uint64_t> &fst_ft, const net_set<uint64_t> &snd_it, const net_set<uint64_t> &snd_ft, bool subtract = false) {
    // float
    auto carry = false;
    auto a_len = fst_ft.length > snd_ft.length ? fst_ft.length : snd_ft.length,
         a_cnt = a_len,
         p_len = a_len;
    auto a_ptr = ptr_init<uint64_t>(p_len);
    auto a_seg = net_ptr_base<uint64_t>{nullptr, 0};
    while (a_len--) {
        auto ans_tmp = dec_add_sub(carry,
                       (fst_ft.length > a_len ? fst_ft[a_len] : 0),
                       (snd_ft.length > a_len ? snd_ft[a_len] : 0),
                       subtract);
        if (a_cnt < p_len || ans_tmp) a_ptr[--a_cnt] = ans_tmp;
    }
    if (a_cnt < p_len) {
        a_seg.init(p_len - a_cnt);
        ptr_copy(a_seg.ptr_base, a_ptr + a_cnt, a_seg.len);
        ans_ft.pointer = std::move(a_seg);
    } else ans_ft.reset();
    // integer
    a_len = fst_it.length > snd_it.length ? fst_it.length : snd_it.length;
    a_cnt = a_len + 1;
    if (p_len < a_cnt) {
        ptr_reset(a_ptr);
        p_len = a_cnt;
        a_ptr = ptr_init<uint64_t>(p_len);
    }
    a_cnt = 0;
    for (auto i = 0ull; i < a_len; ++i) {
        auto ans_tmp = dec_add_sub(carry, 
                       (fst_it.length > i ? fst_it[i] : 0),
                       (snd_it.length > i ? snd_it[i] : 0), 
                       subtract);
        a_ptr[a_cnt++] = ans_tmp;
    }
    if (carry) a_ptr[a_cnt++] = carry;
    else while (a_cnt && !a_ptr[a_cnt - 1]) --a_cnt;
    if (a_cnt) {
        a_seg.init(a_cnt);
        ptr_copy(a_seg.ptr_base, a_ptr, a_cnt);
        ans_it.pointer = std::move(a_seg);
    } else ans_it.reset();
    ptr_reset(a_ptr);
}

// unsafe function
callback_dec_arg net_set<arg> dec_coe(uint64_t &ft_len, const net_set<uint64_t> &it, const net_set<uint64_t> &ft) {
    auto dig_cnt = 0ull,
         seg_cnt = dig_cnt;
    auto seg_ans = ptr_init<arg>(NEUNET_DEC_DIG_MAX * (it.length + ft.length));
    for (auto i = ft.length; i; --i) {
        auto seg_tmp = ft[i - 1];
        while (seg_tmp) {
            auto dig_tmp = seg_tmp % 10;
            seg_tmp     /= 10;
            ++dig_cnt;
            if (seg_cnt || dig_tmp) seg_ans[seg_cnt++] = dig_tmp;
        }
        while (dig_cnt++ < NEUNET_DEC_DIG_MAX) seg_ans[seg_cnt++] = 0;
        dig_cnt = 0;
    }
    ft_len = seg_cnt;
    for (auto i = 0; i < it.length; ++i) {
        auto seg_tmp = it[i];
        while (seg_tmp) {
            seg_ans[seg_cnt++] = seg_tmp % 10;
            seg_tmp /= 10;
            ++dig_cnt;
        }
        if (i + 1 < it.length) {
            while (dig_cnt++ < NEUNET_DEC_DIG_MAX) seg_ans[seg_cnt++] = 0;
            dig_cnt = 0;
        }
    }
    net_ptr_base<arg>ans_ptr;
    ans_ptr.init(seg_cnt);
    ptr_copy(ans_ptr.ptr_base, seg_ans, ans_ptr.len);
    ptr_reset(seg_ans);
    return std::move(ans_ptr);
}

// return factor value
uint64_t dec_fft_coe(net_set<std::complex<long double>> &fst_coe, net_set<std::complex<long double>> &snd_coe) {
    auto ans = fst_coe.length + snd_coe.length - 1,
         len = 1ull << num_bit_cnt(ans - 1);
    fst_coe.init(len);
    snd_coe.init(len);
    return ans;
}

net_set<uint64_t> dec_fft_rev(uint64_t len) {
    net_set<uint64_t> ans(len);
    for (int i = 0; i < len; ++i) ans[i] = (ans[i >> 1] >> 1) | ((i & 1) << (std::max(uint64_t(std::log2(len - 1) + 0.5), 1ull) - 1));
    return ans;
}

void dec_fft(net_set<std::complex<long double>> &coe, const net_set<uint64_t> &rev, bool inv = false) {
    for (auto i = 0ull; i < coe.length; ++i) if (i < rev[i]) std::swap(coe[i], coe[rev[i]]);
    auto half = inv ? -1 : 1;
    for (auto i = 1ull; i < coe.length; i <<= 1) {
        std::complex<long double> w{std::cos(NEUNET_DEC_PI / i),
                                    std::sin(NEUNET_DEC_PI / i) * half};
        for (auto j = 0ull; j < coe.length; j += (i << 1)) {
            std::complex<long double> w_x{1, 0};
            for (auto k = 0ull; k < i; ++k, w_x *= w) {
                auto x = coe[j + k],
                     y = w_x * coe[i + j + k];
                coe[j + k]     = x + y;
                coe[i + j + k] = x - y;
            }
        }
    }
}

// unsigned segment number
void dec_fft_mul(uint64_t &precision, net_set<uint64_t> &ans_it, net_set<uint64_t> &ans_ft, const net_set<uint64_t> &fst_it, const net_set<uint64_t> &fst_ft, const net_set<uint64_t> &snd_it, const net_set<uint64_t> &snd_ft) {
    uint64_t fst_ft_len = 0,
             snd_ft_len = 0;
    auto fst_coe = dec_coe<std::complex<long double>>(fst_ft_len, fst_it, fst_ft),
         snd_coe = dec_coe<std::complex<long double>>(snd_ft_len, snd_it, snd_ft);
    dec_fft_coe(fst_coe, snd_coe);
    auto rev_idx = dec_fft_rev(fst_coe.length);
    dec_fft(fst_coe, rev_idx);
    dec_fft(snd_coe, rev_idx);
    net_set<std::complex<long double>> ans_seg(rev_idx.length);
    net_ptr_base<uint64_t> ans_seg_ptr;
    auto carry  = 0ull,
         dg_cnt = 0ull,
         sg_tmp = 0ull,
         pw_cnt = 1ull,
         ft_dig = fst_ft_len + snd_ft_len;
    for (auto i = 0ull; i < rev_idx.length; ++i) ans_seg[i] = fst_coe[i] * snd_coe[i];
    dec_fft(ans_seg, rev_idx, true);
    dg_cnt    = ft_dig % NEUNET_DEC_DIG_MAX;
    precision = ft_dig;
    if (dg_cnt) {
        dg_cnt  = NEUNET_DEC_DIG_MAX - dg_cnt;
        ft_dig += dg_cnt;
        while (dg_cnt--) pw_cnt *= 10;
        dg_cnt = 0;
    }
    if (ft_dig) dg_cnt = ans_seg.length;
    else ans_ft.reset();
    auto ans_ptr = ptr_init<uint64_t>(rev_idx.length);
    for (auto i = 0ull; i < ans_seg.length; ++i) {
        carry  += (int)(ans_seg[i].real() / ans_seg.length + 0.5);
        sg_tmp += carry % 10 * pw_cnt;
        pw_cnt *= 10;
        carry  /= 10;
        if ((i + 1) == ans_seg.length) while (pw_cnt < NEUNET_DEC_SEG_MAX) {
            if (carry) {
                sg_tmp += carry % 10 * pw_cnt;
                carry  /= 10;
            }
            pw_cnt *= 10;
        }
        if (pw_cnt != NEUNET_DEC_SEG_MAX) continue;
        if (ft_dig) {
            // float
            ft_dig -= NEUNET_DEC_DIG_MAX;
            if (dg_cnt < ans_seg.length || sg_tmp) ans_ptr[--dg_cnt] = sg_tmp;
            if ((i + 1) == ans_seg.length) while (ft_dig) {
                if (carry) ans_ptr[--dg_cnt] = carry;
                else ans_ptr[--dg_cnt] = 0;
                ft_dig -= NEUNET_DEC_DIG_MAX;
            }
            if (!ft_dig) {
                // change to integer
                ans_seg_ptr.init(ans_seg.length - dg_cnt);
                ptr_copy(ans_seg_ptr.ptr_base, ans_ptr + dg_cnt, ans_seg_ptr.len);
                ans_ft.pointer = std::move(ans_seg_ptr);
                dg_cnt             = 0;
            }
        } else ans_ptr[dg_cnt++] = sg_tmp; // integer
        pw_cnt = 1;
        sg_tmp = 0;
    }
    if (carry) ans_ptr[dg_cnt++] = carry;
    else while (dg_cnt && !ans_ptr[dg_cnt - 1]) --dg_cnt;
    if (dg_cnt) {
        ans_seg_ptr.init(dg_cnt);
        ptr_copy(ans_seg_ptr.ptr_base, ans_ptr, dg_cnt);
        ans_it.pointer = std::move(ans_seg_ptr);
    } else ans_it.reset();
    ptr_reset(ans_ptr);
}

void dec_mul(uint64_t &precision, net_set<uint64_t> &ans_it, net_set<uint64_t> &ans_ft, const net_set<uint64_t> &fst_it, const net_set<uint64_t> &fst_ft, const net_set<uint64_t> &snd_it, const net_set<uint64_t> &snd_ft) {
    uint64_t fst_ft_len = 0,
             snd_ft_len = 0;
    auto     fst_coe    = dec_coe<uint8_t>(fst_ft_len, fst_it, fst_ft),
             snd_coe    = dec_coe<uint8_t>(snd_ft_len, snd_it, snd_ft);
    net_set<uint64_t>      ans_seg(fst_coe.length + snd_coe.length - 1);
    net_ptr_base<uint64_t> ans_seg_ptr;
    auto carry  = 0ull,
         dg_cnt = 0ull,
         sg_tmp = 0ull,
         pw_cnt = 1ull,
         ft_dig = fst_ft_len + snd_ft_len;
    for (auto i = 0ull; i < fst_coe.length; ++i) for (auto j = 0ull; j < snd_coe.length; ++j) ans_seg[i + j] += fst_coe[i] * snd_coe[j];
    dg_cnt    = ft_dig % NEUNET_DEC_DIG_MAX;
    precision = ft_dig;
    if (dg_cnt) {
        dg_cnt  = NEUNET_DEC_DIG_MAX - dg_cnt;
        ft_dig += dg_cnt;
        while (dg_cnt--) pw_cnt *= 10;
        dg_cnt = 0;
    }
    if (ft_dig) dg_cnt = ans_seg.length;
    else ans_ft.reset();
    auto ans_ptr = ptr_init<uint64_t>(ans_seg.length);
    for (auto i = 0ull; i < ans_seg.length; ++i) {
        carry  += ans_seg[i];
        sg_tmp += carry % 10 * pw_cnt;
        pw_cnt *= 10;
        carry  /= 10;
        if ((i + 1) == ans_seg.length) while (pw_cnt < NEUNET_DEC_SEG_MAX) {
            if (carry) {
                sg_tmp += carry % 10 * pw_cnt;
                carry  /= 10;
            }
            pw_cnt *= 10;
        }
        if (pw_cnt != NEUNET_DEC_SEG_MAX) continue;
        if (ft_dig) {
            // float
            ft_dig -= NEUNET_DEC_DIG_MAX;
            if (dg_cnt < ans_seg.length || sg_tmp) ans_ptr[--dg_cnt] = sg_tmp;
            if ((i + 1) == ans_seg.length) while (ft_dig) {
                if (carry) ans_ptr[--dg_cnt] = carry;
                else ans_ptr[--dg_cnt] = 0;
                ft_dig -= NEUNET_DEC_DIG_MAX;
            }
            if (!ft_dig) {
                // change to integer
                ans_seg_ptr.init(ans_seg.length - dg_cnt);
                ptr_copy(ans_seg_ptr.ptr_base, ans_ptr + dg_cnt, ans_seg_ptr.len);
                ans_ft.pointer = std::move(ans_seg_ptr);
                dg_cnt             = 0;
            }
        } else ans_ptr[dg_cnt++] = sg_tmp; // integer
        pw_cnt = 1;
        sg_tmp = 0;
    }
    if (carry) ans_ptr[dg_cnt++] = carry;
    else while (dg_cnt && !ans_ptr[dg_cnt - 1]) --dg_cnt;
    if (dg_cnt) {
        ans_seg_ptr.init(dg_cnt);
        ptr_copy(ans_seg_ptr.ptr_base, ans_ptr, dg_cnt);
        ans_it.pointer = std::move(ans_seg_ptr);
    } else ans_it.reset();
    ptr_reset(ans_ptr);
}

int8_t dec_div(net_set<int8_t> &divd, const net_set<uint8_t> &divr) {
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
        if (divd_tmp[divd_tmp.length - 1]) divd.init(divd_tmp.length + 1, false);
        else divd.init(divd_tmp.length + 1, false);
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

void dec_div(net_set<uint64_t> &ans_it, net_set<uint64_t> &ans_ft, const net_set<uint64_t> &divd_it, const net_set<uint64_t> &divd_ft, const net_set<uint64_t> &divr_it, const net_set<uint64_t> &divr_ft, uint64_t precision) {
    auto dd_i_len = 0ull,
         dr_i_len = 0ull,
         it_seg_c = 1ull;
    auto divd_coe = dec_coe<int8_t>(dd_i_len, divd_it, divd_ft);
    auto divr_coe = dec_coe<uint8_t>(dr_i_len, divr_it, divr_ft);
    divd_coe.reverse();
    divr_coe.reverse();
    dd_i_len = divd_coe.length - dd_i_len;
    dr_i_len = divr_coe.length - dr_i_len;
    if (dd_i_len >= dr_i_len) it_seg_c = dd_i_len - dr_i_len + 1;
    net_set<int8_t> ans_set(it_seg_c + precision + 2);
    net_set<uint64_t> ans_seg(it_seg_c + precision);
    auto idx_tool = 0ull;
    if (dd_i_len < dr_i_len) idx_tool = dr_i_len - dd_i_len;
    for (auto i = idx_tool; i < ans_set.length - 2; ++i) ans_set[i] = dec_div(divd_coe, divr_coe);
    auto seg_cnt = it_seg_c % NEUNET_DEC_DIG_MAX;
    it_seg_c    /= NEUNET_DEC_DIG_MAX;
    if (seg_cnt) ++it_seg_c;
    else seg_cnt = NEUNET_DEC_DIG_MAX;
    idx_tool = 0;
    for (auto i = 0ull; i < ans_set.length - 2; ++i) {
        ans_seg[idx_tool] *= 10;
        ans_seg[idx_tool] += ans_set[i];
        if (--seg_cnt) continue;
        auto nidx = i + 1;
        while (ans_set[nidx] < 0) {
            --ans_seg[idx_tool];
            ans_set[nidx] += 10;
        }
        if (!ans_set[nidx]) {
            auto carry = 0;
            while (ans_set[nidx + 1] < 0) {
                ans_set[nidx + 1] += 10;
                --carry;
            }
            while (carry < 0) {
                --ans_seg[idx_tool];
                carry += 10;
            }
            ans_set[nidx] = carry;
        }
        ++idx_tool;
        seg_cnt = NEUNET_DEC_DIG_MAX;
    }
    if (seg_cnt) {
        while (seg_cnt--) ans_seg[idx_tool] *= 10;
        seg_cnt  = idx_tool + 1;
    } else seg_cnt = idx_tool;
    idx_tool = 0;
    while (idx_tool < it_seg_c && !ans_seg[idx_tool]) ++idx_tool;
    while (it_seg_c < seg_cnt && !ans_seg[seg_cnt - 1]) --seg_cnt;
    if (idx_tool < it_seg_c) {
        ans_it = ans_seg.sub_set(idx_tool, it_seg_c - 1);
        ans_it.reverse();
    } else ans_it.reset();
    if (it_seg_c < seg_cnt) ans_ft = ans_seg.sub_set(it_seg_c, seg_cnt - 1);
    else ans_ft.reset();
}

class net_decimal {
public:
    uint64_t precision = 0;
    bool     modulus   = false;

    static uint64_t accuracy;
    static bool     fft_mode;

    __declspec(property(get = dec_ft_dig_cnt,
                        put = dec_truncate)) uint64_t    float_digit_count;
    __declspec(property(get = to_string))    std::string string_format;
    // unsigned format
    __declspec(property(get = to_integer))   uint64_t    integer_format;
    __declspec(property(get = to_float))     long double float_point_format;
    __declspec(property(get = dec_int))      net_decimal integer_part;
    __declspec(property(get = dec_float))    net_decimal float_point_part;
    __declspec(property(get = dec_abs))      net_decimal absolute;
    __declspec(property(get = dec_recip))    net_decimal reciprocal;

protected:
    void value_assign(const net_decimal &src) {
        precision = src.precision;
        modulus   = src.modulus;
    }

    void value_copy(const net_decimal &src) {
        value_assign(src);
        val[0] = src.val[0];
        val[1] = src.val[1];
        sgn    = src.sgn;
    }

    void value_move(net_decimal &&src) {
        value_assign(src);
        val[0] = std::move(src.val[0]);
        val[1] = std::move(src.val[1]);
        sgn    = src.sgn;
        src.reset();
    }

    bool dec_is_zero() const { return !(val[it].length || val[ft].length); }

    bool dec_is_one() const { return !val[ft].length && val[it].length == 1 && val[it][0] == 1; }

    static net_decimal dec_ln_4() {
        net_decimal b = 1,
                    c = 2,
                    s = "0.6",
                    o = "0.36",
                    p = 0.,
                    q = b,
                    m = p,
                    n = q;
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
            n  = p / q;
        }
        return 2 * n;
    }
    net_decimal dec_ln_proto() const {
        net_decimal b = 1,
                    p = 0.,
                    q = b,
                    o = p,
                    u = *this - 1,
                    v = *this + 1,
                    r = u * u,
                    s = v * v,
                    m = q,
                    n = p;
        while (m != n) {
            // std::cout << n << std::endl;
            m = std::move(n);
            if (q == v) p += u;
            else {
                o  = b * v;
                p  = o * p + u * q;
                q *= o;
            }
            u *= r;
            v *= s;
            b += 2;
            n  = p / q;
        }
        return n * 2;
    }
    // ln(x + 1) = (-1 < x <= 1)
    net_decimal dec_ln_2() const {
        net_decimal b = 1,
                    p = 0.,
                    q = b,
                    x = *this - 1,
                    o = x,
                    c = b,
                    m = q,
                    n = p;
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
            n     = p / q;
            ++b;
        }
        return n;
    }

   net_decimal dec_sin_cos_iter(net_decimal &u, net_decimal &b, bool print_iter = false) const {
        net_decimal y = 0.,
                    c = 1.,
                    m = c,
                    v = c,
                    q = c,
                    p = y,
                    n = y,
                    o = (*this) * (*this);
        while (m != n) {
            if (print_iter) std::cout << n << std::endl;
            m = std::move(n);
            if (q == v) p += c * u;
            else {
                p  = v * p + c * q * u;
                q *= v;
            }
            u *= o;
            n  = p / q;
            for (auto i = 0; i < 2; ++i) v *= (++b);
            c.sgn = !c.sgn;
        }
        return n;
    }
    
public:
    net_decimal(const long double src = .0l) { dec_init(precision, sgn, val[it], val[ft], src); }
    net_decimal(const char *src) {
        if (!src) {
            std::cerr << "Null pointer." << std::endl;
            return;
        }
        uint64_t str_len = 0,
                 dot_idx = 0;
        if (!dec_init(precision, sgn, val[it], val[ft], src)) std::cerr << "Invalid initializer." << std::endl;
    }
    net_decimal(const net_decimal &src) { value_copy(src); }
    net_decimal(net_decimal &&src) { value_move(std::move(src)); }

    std::string to_string() const {
        auto ft_len = 0ull,
             cnt    = 0ull;
        auto coe    = dec_coe<uint8_t>(ft_len, val[it], val[ft]);
        auto it_len = coe.length - ft_len;
        std::string ans;
        if (sgn) ans.push_back('-');
        if (ft_len == coe.length) ans.push_back('0');
        for (auto i = coe.length; i; --i) {
            if (cnt++ == it_len) ans.push_back('.');
            ans.push_back(coe[i - 1] + '0');
        }
        return ans;
    }

    uint64_t to_integer() const {
        if (val[it].length == 1) return val[it][0];
        return 0;
    }

    long double to_float() const {
        if (val[it].length == 1 && val[ft].length == 1) {
            auto tmp = val[ft][0];
            auto ans = (long double)to_integer(),
                 flt = 0.l;
            while (!(tmp % 10)) tmp /= 10;
            while (tmp) {
                flt += tmp % 10;
                flt /= 10;
                tmp /= 10;
            }
            ans += flt;
            if (sgn) return (-1) * ans;
            else return ans;
        }
        return 0;
    }

    void reset() {
        val[0].reset();
        val[1].reset();
        sgn = false;
    }

    virtual ~net_decimal() { reset(); }

protected:
    /* 1|1010|0100|0131|0800|1920|0001.0001|0400|0021|7080|0400|005
     * {1,1920,800,131,100,1010,1}     {1,400,21,7080,400,50}
     */
    net_set<uint64_t> val[2];

    static constexpr bool it = false,
                          ft = true;

    bool sgn = false;

public:
    static net_decimal dec_pi() {
        net_decimal p = 0.,
                    q = 1,
                    a = 0.25, // +2
                    b = 2,    // +4
                    c = 5,    // +8
                    d = 6,    // +8
                    o = 0.0625,
                    s = q,
                    m = p,
                    n = q;
        while (m != n) {
            // std::cout << n << std::endl;
            m = std::move(n);
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
            n = p / q;
        }
        return n;
    }

    net_decimal dec_abs() const {
        auto ans = *this;
        if (ans.sgn) ans.sgn = false;
        return ans;
    }

    uint64_t dec_ft_dig_cnt() const {
        if (!val[ft].length) return 0;
        auto ans = val[ft].length * NEUNET_DEC_DIG_MAX;
        auto seg = val[ft][val[ft].length - 1];
        while (!(seg % 10)) {
            seg /= 10;
            --ans;
        }
        return ans;
    }

    void dec_truncate(uint64_t acc) {
        auto dig_cnt = dec_ft_dig_cnt();
        if (acc < dig_cnt) acc += 2;
        else return;
        net_decimal carry = 1;
        for (auto i = 0ull; i < acc - 1; ++i) carry *= 0.1;
        uint64_t tgt_idx = acc / NEUNET_DEC_DIG_MAX,
                 seg_len = acc % NEUNET_DEC_DIG_MAX,
                 dig_idx = NEUNET_DEC_DIG_MAX;
        if (seg_len) ++tgt_idx;
        else seg_len = dig_idx;
        val[ft].init(tgt_idx--);
        auto curr_seg = val[ft][tgt_idx],
             seg_pow  = 1ull,
             last_dig = 0ull;
        // get last 2 digits after accuracy digits
        while (dig_idx-- > seg_len) {
            last_dig          = seg_pow * (curr_seg % 10);
            curr_seg         /= 10;
            seg_pow          *= 10;
            val[ft][tgt_idx] -= last_dig;
        }
        last_dig = curr_seg % 10;
        curr_seg /= 10;
        // last digit of the 2;
        if (last_dig >= 5) {
            // carry
            *this += carry;
            ++curr_seg;
        }
        last_dig *= seg_pow;
        seg_pow  *= 10;
        if (dig_idx) {
            // not the end of current segment
            val[ft][tgt_idx] -= last_dig;
            last_dig          = curr_seg % 10;
            if (last_dig < 5) val[ft][tgt_idx] -= last_dig * seg_pow;
            if (!val[ft][tgt_idx]) val[ft].init(tgt_idx);
        } else {
            // current segment ends
            val[ft].init(tgt_idx--);
            last_dig = val[ft][tgt_idx] % 10;
            if (last_dig < 5) val[ft][tgt_idx] -= last_dig;
        }
        // mark
        precision = acc - 2;
    }

    net_decimal dec_int() const {
        net_decimal ans;
        ans.sgn     = sgn;
        ans.val[it] = val[it];
        return ans;
    }

    net_decimal dec_float() const {
        net_decimal ans;
        ans.sgn     = sgn;
        ans.val[ft] = val[ft];
        return ans;
    }

    net_decimal dec_recip() const {
        if (dec_is_zero()) return 0.;
        if (dec_is_one()) return sgn ? (-1) : 1;
        auto trunc_cnt   = accuracy ? accuracy : precision;
        net_decimal curr = .1, prev = .0, coe = 2;
        while (curr * (*this) > 1) curr *= 0.1;
        // auto ans = curr;
        while (curr != prev) {
            prev = std::move(curr);
            curr = prev * (coe - prev * (*this));
            // curr = ans;
            curr.dec_truncate(trunc_cnt);
            // std::cout << curr << std::endl;
        }
        return curr;
    }

    net_decimal dec_ln() const {
        if (dec_is_zero() || dec_is_one() || sgn) return 0.;
        net_decimal ln_4 = 0.,
                    base = *this;
        while (base > 1) {
            base *= 0.25;
            ++ln_4;
        }
        ln_4 *= dec_ln_4();
        return base.dec_ln_2() + ln_4;
    }

    net_decimal dec_exp() const {
        net_decimal c = 1,
                    p = c,
                    q = c,
                    v = c,
                    u = *this,
                    m = 0.,
                    n = 1;
        while (m != n) {
            m = std::move(n);
            if (q == v) p += u;
            else {
                p  = v * p + u * q;
                q *= v;
            }
            ++c;
            u *= *this;
            v *= c;
            n  = p / q;
        }
        return n;
    }

    net_decimal dec_sin(bool print_iter = false) const {
        net_decimal u = (*this),
                    b = 1;
        return dec_sin_cos_iter(u, b, print_iter);
    }

    net_decimal dec_cos(bool print_iter = false) const {
        net_decimal u = 1,
                    b = 0.;
        return dec_sin_cos_iter(u, b, print_iter);
    }

    net_decimal dec_pow(const net_decimal &times) const {
        if (dec_is_zero()) return 0.;
        if ((!sgn && dec_is_one()) || times.dec_is_zero()) return 1;
        if (times.dec_is_one()) {
            if (times.sgn) return 1 / (*this);
            return *this;
        }
        if (times.val[ft].length) if (sgn) {
            net_decimal ans = 0.,
                        cnt = 1.,
                        pia = dec_pi() * times,
                        prd = ans;
            net_decimal res[NEUNET_CACHE_LEN];
            uint64_t    res_cnt = 0;
            do {
                prd = cnt * pia;
                ans = prd.dec_sin();
                ans.dec_truncate(precision - 1);
                // std::cout << ans << '\n' << std::endl;
                if (ans.dec_is_zero()) {
                    // std::cout << std::endl;
                    prd = prd.dec_cos();
                    // std::cout << prd << std::endl;
                    // prd.dec_truncate(times.accuracy - 2);
                    // std::cout << prd << std::endl;
                    return dec_abs().dec_pow(times) * prd;
                }
                res[res_cnt++] = (ans);
                cnt += 2;
            } while (res_cnt == 1 || ans != res[0]);
            return ans;
        } else return (times * dec_ln()).dec_exp();
        else {
            auto ans = *this,
                 bas = ans,
                 cnt = times,
                 pwr = net_decimal(2),
                 dcy = pwr;
            cnt.sgn  = false;
            while (dec_compare(pwr.val[it], pwr.val[ft], times.val[it], times.val[ft]) == NEUNET_DEC_CMP_LES) {
                ans *= ans;
                pwr *= dcy;
            }
            cnt -= pwr * 0.5;
            while (!cnt.dec_is_zero()) {
                ans *= bas;
                --cnt;
            }
            if (times.sgn) return 1 / ans;
            return ans;
        }
    }

    net_decimal &operator=(const net_decimal &src) {
        value_copy(src);
        return *this;
    }
    net_decimal &operator=(net_decimal &&src) {
        value_move(std::move(src));
        return *this;
    }

    friend bool operator==(const net_decimal &fst, const net_decimal &snd) {
        if (&fst == &snd) return true;
        return fst.sgn == snd.sgn && dec_compare(fst.val[it], fst.val[ft], snd.val[it], snd.val[ft]) == NEUNET_DEC_CMP_EQL;
    }
    friend bool operator!=(const net_decimal &fst, const net_decimal &snd) { return !(fst == snd); }

    friend bool operator>(const net_decimal &fst, const net_decimal &snd) {
        if (&fst == &snd) return false;
        if (fst.sgn != snd.sgn)
            if (snd.sgn) return true;
            else return false;
        else {
            auto dec_cmp = dec_compare(fst.val[it], fst.val[ft], snd.val[it], snd.val[ft]);
            if (fst.sgn && dec_cmp == NEUNET_DEC_CMP_LES || !fst.sgn && dec_cmp == NEUNET_DEC_CMP_GTR) return true;
            else return false;
        }
    }
    friend bool operator>=(const net_decimal &fst, const net_decimal &snd) { return fst > snd || fst == snd; }

    friend bool operator<(const net_decimal &fst, const net_decimal &snd) { return !(fst >= snd); }
    friend bool operator<=(const net_decimal &fst, const net_decimal &snd) { return !(fst > snd); }

    net_decimal operator+() { return *this; }
    friend net_decimal operator+(const net_decimal &fst, const net_decimal &snd) {
        if (fst.dec_is_zero()) return snd;
        if (snd.dec_is_zero()) return fst;
        net_decimal ans;
        if (fst.sgn == snd.sgn) {
            dec_add_sub(ans.val[it], ans.val[ft], fst.val[it], fst.val[ft], snd.val[it], snd.val[ft]);
            if (fst.sgn) ans.sgn = true;
        } else {
            auto dec_cmp = dec_compare(snd.val[it], snd.val[ft], fst.val[it], fst.val[ft]);
            if (dec_cmp == NEUNET_DEC_CMP_GTR) {
                dec_add_sub(ans.val[it], ans.val[ft], snd.val[it], snd.val[ft], fst.val[it], fst.val[ft], true);
                if (snd.sgn) ans.sgn = true;
            }
            if (dec_cmp == NEUNET_DEC_CMP_LES) {
                dec_add_sub(ans.val[it], ans.val[ft], fst.val[it], fst.val[ft], snd.val[it], snd.val[ft], true);
                if (fst.sgn) ans.sgn = true;
            }
        }
        ans.precision = dec_add_prec(fst.precision, snd.precision);
        return ans;
    }
    callback_dec_arg friend void operator+=(arg &fst, const net_decimal &snd) { fst += snd.to_float(); }
    void operator+=(const net_decimal &src) { *this = *this + src; }
    net_decimal &operator++() {
        *this += 1;
        return *this;
    }
    net_decimal operator++(int) {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    net_decimal operator-() { return 0. - *this; }
    friend net_decimal operator-(const net_decimal &fst, const net_decimal &snd) {
        if (snd.dec_is_zero()) return fst;
        if (fst.dec_is_zero()) {
            auto tmp = snd;
            tmp.sgn  = !snd.sgn;
            return tmp;
        }
        net_decimal ans;
        if (fst.sgn == snd.sgn) {
            auto dec_cmp = dec_compare(fst.val[it], fst.val[ft], snd.val[it], snd.val[ft]);
            if (dec_cmp == NEUNET_DEC_CMP_GTR) {
                dec_add_sub(ans.val[it], ans.val[ft], fst.val[it], fst.val[ft], snd.val[it], snd.val[ft], true);
                if (fst.sgn) ans.sgn = true;
            }
            if (dec_cmp == NEUNET_DEC_CMP_LES) {
                dec_add_sub(ans.val[it], ans.val[ft], snd.val[it], snd.val[ft], fst.val[it], fst.val[ft], true);
                if (!fst.sgn) ans.sgn = true;
            }
        } else {
            dec_add_sub(ans.val[it], ans.val[ft], fst.val[it], fst.val[ft], snd.val[it], snd.val[ft]);
            if (fst.sgn) ans.sgn = true;
        }
        ans.precision = dec_add_prec(fst.precision, snd.precision);
        return ans;
    }
    callback_dec_arg friend void operator-=(arg &fst, const net_decimal &snd) { fst += snd.to_float(); }
    void operator-=(const net_decimal &src) { *this = *this - src; }
    net_decimal &operator--() {
        *this -= 1;
        return *this;
    }
    net_decimal operator--(int) {
        auto tmp = *this;
        --(*this);
        return tmp;
    }

    friend net_decimal operator*(const net_decimal &fst, const net_decimal &snd) {
        net_decimal ans;
        if (fst.dec_is_zero() || snd.dec_is_zero()) return ans;
        if (fst.dec_is_one()) {
            ans = snd;
            if (fst.sgn) ans.sgn = !ans.sgn;
            return ans;
        }
        if (snd.dec_is_one()) {
            ans = fst;
            if (snd.sgn) ans.sgn = !ans.sgn;
            return ans;
        }
        if (fft_mode) dec_fft_mul(ans.precision, ans.val[it], ans.val[ft], fst.val[it], fst.val[ft], snd.val[it], snd.val[ft]);
        else dec_mul(ans.precision, ans.val[it], ans.val[ft], fst.val[it], fst.val[ft], snd.val[it], snd.val[ft]);
        ans.sgn = fst.sgn != snd.sgn;
        return ans;
    }
    callback_dec_arg friend void operator*=(arg &fst, const net_decimal &snd) { fst *= snd.to_float(); }
    void operator*=(const net_decimal &src) { *this = *this * src; }

    friend net_decimal operator/(const net_decimal &fst, const net_decimal &snd) {
        assert(!snd.dec_is_zero());
        net_decimal ans;
        if (fst.dec_is_zero()) return ans;
        if (snd.dec_is_one()) {
            if (fst.sgn == snd.sgn) return fst;
            ans     = fst;
            ans.sgn = true;
            return ans;
        }
        ans.sgn = snd.sgn != fst.sgn;
        dec_div(ans.val[it], ans.val[ft], fst.val[it], fst.val[ft], snd.val[it], snd.val[ft], accuracy + 1);
        ans.precision = accuracy;
        return ans;
    }
    callback_dec_arg friend void operator/=(arg &fst, const net_decimal &snd) { fst /= snd.to_float(); }
    void operator/=(const net_decimal &src) { *this = *this / src; }

    friend net_decimal operator%(const net_decimal &fst, const net_decimal &snd) {
        if (fst.val[ft].length || snd.val[ft].length) return 0.;
        auto ans = (fst / snd).dec_int();
        if (snd.modulus && ans.sgn) --ans;
        return fst - snd * ans;
    }
    void operator%=(const net_decimal &src) { *this = *this % src; }

    friend std::ostream &operator<<(std::ostream &out, const net_decimal &src) { return out << src.to_string(); }
};
uint64_t net_decimal::accuracy = 32;
bool     net_decimal::fft_mode = false;

NEUNET_END

_STD_BEGIN

neunet::net_decimal abs(const neunet::net_decimal &src) { return src.absolute; }

_STD_END
