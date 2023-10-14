NEUNET_BEGIN

class net_decimal {
public:
    bool modulus_mode = false;

    uint64_t division_precision = 32,
             binary_bit_count   = 0;

protected:
    void value_assign(const net_decimal &src) {
        sgn                = src.sgn;
        modulus_mode       = src.modulus_mode;
        division_precision = src.division_precision;
    }

    void value_copy(const net_decimal &src) {
        value_assign(src);
        num = src.num;
        den = src.den;
    }

    void value_move(net_decimal &&src) {
        value_assign(src);
        num = std::move(src.num);
        den = std::move(src.den);
    }

    bool is_zero() const { return dec_is_zero(num); }

    bool is_one() const { return dec_comp(den, num) == NEUNET_DEC_CMP_EQL; }

    bool is_dec() const { return dec_is_zero(den); }

    bool is_integer() {
        auto zero_flag = is_zero(),
             one_flag  = is_one();
        if (zero_flag || one_flag) {
            den.reset();
            if (one_flag) num = dec_init(sgn, 1);
            return true;
        }
        if (is_dec()) return !num.ft.length;
        reduct();
        return dec_is_one(den);
    }

    void frac_norm() {
        if (is_dec()) den = dec_init(red, 1);
        net_decimal_data e10;
        if (num.ft.length) {
            num = dec_e10(e10, num);
            den = dec_mul(e10, den);
        }
        if (den.ft.length) {
            den = dec_e10(e10, den);
            num = dec_mul(e10, num);
        }
    }

    bool sin_cos_period(net_decimal &src) const {
        net_decimal pi_2;
        pi_2.num = dec_init(pi_2.sgn, 2);
        pi_2.num = dec_mul(pi_2.num, net_decimal_pi.get(division_precision * 2));
        if (*this <= pi_2) return false;
        auto tmp {*this / pi_2};
        dec_float_part(tmp.num, dec_div(tmp.num, tmp.den, 4));
        tmp.den.reset();
        src = *this - tmp * pi_2;
        src.division_precision = division_precision;
        return true;
    }

public:
    net_decimal() {}
    net_decimal(const net_decimal &src) { value_copy(src); }
    net_decimal(net_decimal &&src) { value_move(std::move(src)); }
    callback_arg net_decimal(const arg &src) {
        static_assert(neunet_dec_init_expr);
        num = dec_init(sgn, src);
    }

    void reduct() {
        if (red) return;
        else red = true;
        if (is_dec() && !num.ft.length) return;
        auto zero_valid = is_zero(),
             one_valid  = is_one();
        if (zero_valid || one_valid) {
            if (one_valid) num = den;
            return;
        }
        if (num.ft.length || den.ft.length) frac_norm();
        auto fst = num,
             snd = den;
        auto res = dec_gcd(fst, snd);
        net_decimal_data fac;
        if (res) fac = std::move(snd);
        else fac = std::move(fst);
        num = dec_rem(num, fac);
        den = dec_rem(den, fac);
    }

    int64_t to_integer() const {
        if (is_dec()) return dec2i(sgn, num);
        net_decimal_data ans_num = num;
        ans_num = dec_rem(ans_num, den);
        return dec2i(sgn, ans_num);
    }

    long double to_float() const {
        if (is_dec()) return dec2f(sgn, num);
        return dec2f(sgn, dec_div(num, den, NEUNET_DEC_VLD_DIG - 2));
    }

    net_decimal float_part(net_decimal &integer_part) const {
        integer_part.reset();
        integer_part.sgn = sgn;
        net_decimal ans;
        ans.sgn = sgn;
        if (is_dec()) {
            ans.num = dec_float_part(integer_part.num, num);
            return ans;
        }
        ans = *this;
        ans.reduct();
        integer_part.num = dec_rem(ans.num, ans.den);
        return ans;
    }

    net_decimal remainder(net_decimal &divr) {
        net_decimal ans = *this;
        dec_rem(ans.num, divr.num);
        if (divr.modulus_mode && sgn != divr.sgn) ans.num = dec_add(ans.sgn, ans.num, ans.sgn, divr.num, divr.sgn);
        return ans;
    }

    callback_dec_arg static void remainder(arg &divd, net_decimal &divr) {
        static_assert(std::is_integral_v<arg>, "Dividend should be an integer.");
        auto divr_arith = divr.to_integer();
        auto divd_sgn   = divd < 0;
        divd = divd % divr_arith;
        if (divr.modulus_mode && divd_sgn != divr.sgn) divd += divr_arith;
        if (divd_sgn == divr.sgn && divd_sgn) divd *= (-1);
    }

    int compare(const net_decimal &src) const {
        if (sgn != src.sgn) return sgn ? NEUNET_DEC_CMP_LES : NEUNET_DEC_CMP_GTR;
        auto this_dec = is_dec(),
             src_dec  = src.is_dec();
        if (this_dec && src_dec) return dec_comp(sgn, dec_comp(num, src.num));
        if (this_dec && !src_dec) return dec_comp(sgn, dec_comp(dec_mul(num, src.den), src.num));
        if (!this_dec && src_dec) return dec_comp(sgn, dec_comp(num, dec_mul(den, src.num)));
        return dec_comp(sgn, dec_comp(dec_mul(num, src.den), dec_mul(src.num, den)));
    }

    void reset() {
        num.reset();
        den.reset();        
        sgn                = false;
        modulus_mode       = false;
        division_precision = 32;
    }

    ~net_decimal() { reset(); }

    net_decimal abs() const {
        if (!sgn) return *this;
        net_decimal ans {*this};
        ans.sgn = false;
        return ans;
    }

    net_decimal ln() const {
        if (is_zero() || is_one() || sgn) return {};
        net_decimal ans;
        auto minu_sgn = false;
        auto minu_ans = dec_ln(minu_sgn, num, division_precision);
        if (is_dec()) {
            ans.num = std::move(minu_ans);
            ans.sgn = minu_sgn;
            return ans;
        }
        auto subt_sgn = false;
        auto subt_ans = dec_ln(subt_sgn, den, division_precision);
        ans.num = dec_sub(ans.sgn, minu_ans, minu_sgn, subt_ans, subt_sgn);
        return ans;
    }

    net_decimal exp() const {
        // TODO optimize procedure
        net_decimal ans;
        if (is_zero()) {
            ans.num = dec_init(ans.sgn, 1);
            return ans;
        }
        ans.num = dec_exp(ans.sgn, num, den, sgn, division_precision);
        return ans;
    }

    net_decimal sin() const {
        // TODO inducing formula
        net_decimal src;
        if (sin_cos_period(src)) return src.sin();
        net_decimal ans;
        auto b  = dec_init(ans.sgn, 1),
             u  = num,
             v  = den;
        ans.num = dec_sin_cos(ans.sgn, num, den, sgn, u, v, b, division_precision);
        return ans;
    }

    net_decimal cos() const {
        // TODO inducing formula
        net_decimal src;
        if (sin_cos_period(src)) return src.cos();
        net_decimal ans;
        auto b = dec_init(ans.sgn, 0),
             u = dec_init(ans.sgn, 1),
             v = u;
        ans.num = dec_sin_cos(ans.sgn, num, den, sgn, u, v, b, division_precision);
        return ans;
    }

    net_decimal pow(const net_decimal &times) const {
        if (is_zero()) return 0.;
        if ((!sgn && is_one()) || times.is_zero()) return 1;
        if (times.is_one()) {
            if (times.sgn) return 1 / (*this);
            return *this;
        }
        if (!sgn) return (times * ln()).exp();
        net_decimal piv {0};
        piv.num = net_decimal_pi.get(division_precision + 1);
        net_decimal ans {0},
                    pia {piv * times},
                    prd {ans};
        net_decimal res[NEUNET_CACHE_LEN];
        uint64_t    res_cnt {0},
                    cnt {1};
        do {
            prd = cnt * pia;
            ans = prd.sin();
            if (ans.is_zero()) {
                prd = prd.cos();
                dec_truncate(prd.num, division_precision);
                if (sgn) return abs().pow(times) * prd;
                return pow(times) * prd;
            }
            res[res_cnt++] = (ans);
            cnt += 2;
        } while (res_cnt == 1 || ans != res[0]);
        return ans;
    }

    void truncate() {
        if (!is_dec()) return;
        dec_truncate(num, division_precision);
    }

protected:
    // numerator / denominator
    net_decimal_data num, 
                     den;

    bool sgn = false,
         red = false;

public:
    net_decimal &operator=(const net_decimal &src) {
        value_copy(src);
        return *this;
    }
    net_decimal &operator=(net_decimal &&src) {
        value_move(std::move(src));
        return *this;
    }

    friend bool operator==(const net_decimal &fst, const net_decimal &snd) { return fst.compare(snd) == NEUNET_DEC_CMP_EQL; }
    friend std::strong_ordering operator<=>(const net_decimal &fst, const net_decimal &snd) {
        auto cmp = fst.compare(snd);
        if (cmp == NEUNET_DEC_CMP_LES) return std::strong_ordering::less;
        if (cmp == NEUNET_DEC_CMP_GTR) return std::strong_ordering::greater;
        return std::strong_ordering::equal;
    }

    net_decimal operator+() { return *this; }
    friend net_decimal operator+(const net_decimal &fst, const net_decimal &snd) {
        net_decimal ans;
        auto fst_dec = fst.is_dec(),
             snd_dec = snd.is_dec();
        if (fst_dec && snd_dec) {
            ans.num = dec_add(ans.sgn, fst.num, fst.sgn, snd.num, snd.sgn);
            return ans;
        }
        if (fst_dec && !snd_dec) {
            ans.den = snd.den;
            ans.num = dec_add(ans.sgn, dec_mul(fst.num, snd.den), fst.sgn, snd.num, snd.sgn);
            return ans;
        }
        if (!fst_dec && snd_dec) {
            ans.den = fst.den;
            ans.num = dec_add(ans.sgn, fst.num, fst.sgn, dec_mul(snd.num, fst.den), snd.sgn);
            return ans;
        }
        ans.den = dec_mul(fst.den, snd.den);
        ans.num = dec_add(ans.sgn, dec_mul(fst.num, snd.den), fst.sgn, dec_mul(snd.num, fst.den), snd.sgn);
        return ans;
    }
    void operator+=(const net_decimal &snd) { *this = *this + snd; }
    callback_dec_arg friend void operator+=(arg &fst, const net_decimal &snd) { fst += snd.to_float(); }
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
    friend net_decimal operator-(const net_decimal &minu, const net_decimal &subt) {
        net_decimal ans;
        auto minu_dec = minu.is_dec(),
             subt_dec = subt.is_dec();
        if (minu_dec && subt_dec) {
            ans.num = dec_sub(ans.sgn, minu.num, minu.sgn, subt.num, subt.sgn);
            return ans;
        }
        if (minu_dec && !subt_dec) {
            ans.den = subt.den;
            ans.num = dec_sub(ans.sgn, dec_mul(minu.num, subt.den), minu.sgn, subt.num, subt.sgn);
            return ans;
        }
        if (!minu_dec && subt_dec) {
            ans.den = minu.den;
            ans.num = dec_sub(ans.sgn, minu.num, minu.sgn, dec_mul(subt.num, minu.den), subt.sgn);
            return ans;
        }
        ans.den = dec_mul(minu.den, subt.den);
        ans.num = dec_sub(ans.sgn, dec_mul(minu.num, subt.den), minu.sgn, dec_mul(subt.num, minu.den), subt.sgn);
        return ans;
    }
    void operator-=(const net_decimal &subt) { *this = *this - subt; }
    callback_dec_arg friend void operator-=(arg &minu, const net_decimal &subt) { minu -= subt.to_float(); }
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
        auto fst_dec = fst.is_dec(),
             snd_dec = snd.is_dec();
        ans.num      = dec_mul(ans.sgn, fst.num, fst.sgn, snd.num, snd.sgn);
        if (fst_dec && snd_dec) return ans;
        if (fst_dec || snd_dec) 
            if (fst_dec) ans.den = snd.den;
            else ans.den = fst.den;
        else ans.den = dec_mul(fst.den, snd.den);
        return ans;
    }
    void operator*=(const net_decimal &snd) { *this = *this * snd; }
    callback_dec_arg friend void operator*=(arg &fst, const net_decimal &snd) { fst *= snd.to_float(); }

    friend net_decimal operator/(const net_decimal &divd, const net_decimal &divr) {
        net_decimal ans;
        auto divd_dec = divd.is_dec(),
             divr_dec = divr.is_dec();
        ans.sgn = divd.sgn != divr.sgn;
        if (divd_dec && divr_dec) {
            ans.num = divd.num;
            ans.den = divr.num;
            return ans;
        }
        if (divd_dec && !divr_dec) {
            ans.den = divr.num;
            ans.num = dec_mul(divd.num, divr.den);
            return ans;
        }
        if (!divd_dec && divr_dec) {
            ans.num = divd.num;
            ans.den = dec_mul(divd.den, divr.num);
            return ans;
        }
        ans.num = dec_mul(divd.num, divr.den);
        ans.den = dec_mul(divd.den, divr.num);
        return ans;
    }
    void operator/=(const net_decimal &divr) { *this = *this / divr; }
    callback_dec_arg friend void operator/=(arg &divd, const net_decimal &divr) { divd *= divr.to_float(); }

    friend net_decimal operator%(const net_decimal &divd, const net_decimal &divr) {
        auto divd_tmp = divd,
             divr_tmp = divr;
        return divd_tmp.remainder(divr_tmp);
    }
    void operator%=(const net_decimal &divr) { *this = *this % divr; }
    callback_dec_arg friend void operator%=(arg &divd, const net_decimal &divr) {
        static_assert(std::is_integral_v<arg>, "Dividend should be an integer.");
        net_decimal divr_tmp = divr;
        remainder(divd, divr_tmp);
    }

    friend std::ostream &operator<<(std::ostream &os, const net_decimal &src) {
        if (src.sgn && !src.is_zero()) os << '-';
        if (src.is_dec()) {
            os << src.num;
            return os;
        }
        auto tmp = dec_div(src.num, src.den, src.division_precision);
        os << tmp;
        return os;
    }
};

net_decimal operator""_d(const char *src, uint64_t len) { return net_decimal(std::string(src)); }
net_decimal operator""_d(long double src) { return net_decimal(src); }
net_decimal operator""_d(uint64_t src) { return net_decimal(src); }

NEUNET_END

_STD_BEGIN

neunet::net_decimal abs(const neunet::net_decimal &_Xx) { return _Xx.abs(); }

neunet::net_decimal pow(const neunet::net_decimal &_Xx, const neunet::net_decimal &_Yx) { return _Xx.pow( _Yx); }

neunet::net_decimal exp(const neunet::net_decimal &_Xx) { return _Xx.exp(); }

neunet::net_decimal log(const neunet::net_decimal &_Xx) { return _Xx.ln(); }

neunet::net_decimal sin(const neunet::net_decimal &_Xx) { return _Xx.sin(); }

neunet::net_decimal cos(const neunet::net_decimal &_Xx) { return _Xx.cos(); }

_STD_END
