NEUNET_BEGIN

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
        return os;
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

std::string dec_to_string(bool sgn, net_decimal_frac &src, uint64_t prec) {
    if (dec_is_one(src.den)) return dec_to_string(sgn, src.num);
    dec_reduct(src);
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

NEUNET_END