NEUNET_BEGIN

class net_decimal {
public:
    bool modulus_mode = false;
    
    uint64_t binary_bit_set     = 0,
             division_precision = 64;

    __declspec(property(get = abs))      net_decimal absolute;
    __declspec(property(get = to_float)) long double float_point_format;
    __declspec(property(get = to_int))   int64_t     integer_format;

protected:
    void value_assign(const net_decimal &src) {
        sgn                = src.sgn;
        red                = src.red;
        modulus_mode       = src.modulus_mode;
        binary_bit_set     = src.binary_bit_set;
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

    template <uint8_t opt_idx> static net_decimal binary_operator(const net_decimal &fst, const net_decimal &snd) {
        net_decimal ans, fst_red, snd_red;
        if (!dec_is_zero(fst.den)) {
            fst_red = fst;
            fst_red.reduct();
            if (!dec_is_zero(fst_red.den)) return ans;
            return binary_operator<opt_idx>(fst_red, snd);
        }
        if (!dec_is_zero(snd.den)) {
            snd_red = fst;
            snd_red.reduct();
            if (!dec_is_zero(snd_red.den)) return ans;
            return binary_operator<opt_idx>(fst, snd_red);
        }
        if (fst.num.ft.length || snd.num.ft.length) return ans;
        auto bit_set = 0ull;
        if (fst.binary_bit_set && snd.binary_bit_set) bit_set = (std::max)(fst.binary_bit_set, snd.binary_bit_set);
        if constexpr (opt_idx == NEUNET_DEC_BIN_OR) ans.num = dec_bit_or(fst.num, snd.num, fst.sgn, snd.sgn, bit_set);
        else if constexpr (opt_idx == NEUNET_DEC_BIN_AND) ans.num = dec_bit_and(fst.num, snd.num, fst.sgn, snd.sgn, bit_set);
        else if constexpr (opt_idx == NEUNET_DEC_BIN_XOR) ans.num = dec_bit_xor(fst.num, snd.num, fst.sgn, snd.sgn, bit_set);
        return ans;
    }

public:
    net_decimal() { num = dec_init(sgn, 0); }
    net_decimal(const net_decimal &src) { value_copy(src); }
    net_decimal(net_decimal &&src) { value_move(std::move(src)); }
    callback_arg net_decimal(const arg &src) {
        static_assert(neunet_dec_init_expr);
        num = dec_init(sgn, src);
    }

    void reduct() {
        if (red) return;
        red = true;
        dec_frac_red(num, den);
    }

    int64_t to_int() { return dec_frac2i(sgn, num, den); }

    long double to_float() { return dec_frac2f(sgn, num, den); }

    net_decimal int_part(net_decimal &float_part) const {
        float_part.num = num;
        float_part.den = den;
        net_decimal ans;
        ans.num = dec_frac_int_part(float_part.num, float_part.den);
        return ans;
    }

    net_decimal abs() const {
        if (!sgn) return *this;
        auto ans = *this;
        ans.sgn  = false;
        return ans;
    }

    void reset() {
        sgn = false;
        red = true;
        num.reset();
        den.reset();
        modulus_mode       = false;
        binary_bit_set     = 0;
        division_precision = 64;
    }

    ~net_decimal() { reset(); }

protected:
    // numerator / denominator
    net_decimal_data num, den;

    bool sgn = false,
         red = true;

public:
    callback_dec_arg explicit operator arg() const {
        if (dec_is_zero(den)) {
            if constexpr (std::is_unsigned_v<neunet_dec_type(arg)>) {
                if (sgn) std::abort();
                return num.it[0];
            }
            if constexpr (std::is_integral_v<neunet_dec_type(arg)>) return dec2i(sgn, num);
            return dec2f(sgn, num);
        }
        net_decimal tmp;
        tmp.num = dec_div(num, den, division_precision);
        return arg(tmp);
    }

    net_decimal &operator=(const net_decimal &src) {
        value_copy(src);
        return *this;
    }
    net_decimal &operator=(net_decimal &&src) {
        value_move(std::move(src));
        return *this;
    }

    friend bool operator==(const net_decimal &fst, const net_decimal &snd) { return dec_frac_comp(fst.num, fst.den, fst.sgn, snd.num, snd.den, snd.sgn) == NEUNET_DEC_CMP_EQL; }
    friend std::strong_ordering operator<=>(const net_decimal &fst, const net_decimal &snd) {
        auto cmp = dec_frac_comp(fst.num, fst.den, fst.sgn, snd.num, snd.den, snd.sgn);
        if (cmp == NEUNET_DEC_CMP_LES) return std::strong_ordering::less;
        if (cmp == NEUNET_DEC_CMP_GTR) return std::strong_ordering::greater;
        return std::strong_ordering::equal;
    }

    net_decimal operator+() { return *this; }
    friend net_decimal operator+(const net_decimal &fst, const net_decimal &snd) {
        net_decimal ans;
        ans.sgn = dec_frac_add(ans.num, ans.den, fst.num, fst.den, fst.sgn, snd.num, snd.den, snd.sgn);
        ans.red = dec_frac1(ans.num, ans.den) || dec_frac0(ans.num, ans.den) || dec_is_zero(ans.den);
        return ans;
    }
    net_decimal &operator+=(const net_decimal &snd) {
        *this = *this + snd;
        return *this;
    }
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
    friend net_decimal operator-(const net_decimal &fst, const net_decimal &snd) {
        net_decimal ans;
        ans.sgn = dec_frac_sub(ans.num, ans.den, fst.num, fst.den, fst.sgn, snd.num, snd.den, snd.sgn);
        ans.red = dec_frac1(ans.num, ans.den) || dec_frac0(ans.num, ans.den) || dec_is_zero(ans.den);
        return ans;
    }
    net_decimal &operator-=(const net_decimal &snd) {
        *this = *this - snd;
        return *this;
    }
    callback_dec_arg friend arg operator-=(arg &fst, const net_decimal &snd) { return fst -= snd.to_float(); }
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
        ans.sgn = dec_frac_mul(ans.num, ans.den, fst.num, fst.den, fst.sgn, snd.num, snd.den, snd.sgn);
        ans.red = dec_frac1(ans.num, ans.den) || dec_frac0(ans.num, ans.den) || dec_is_zero(ans.den);
        return ans;
    }
    net_decimal &operator*=(const net_decimal &snd) {
        *this = *this * snd;
        return *this;
    }
    callback_dec_arg friend arg operator*=(arg &fst, const net_decimal &snd) { return fst *= snd.to_float(); }

    friend net_decimal operator/(const net_decimal &divd, const net_decimal &divr) {
        net_decimal ans;
        ans.sgn = dec_frac_div(ans.num, ans.den, divd.num, divd.den, divd.sgn, divr.num, divr.den, divr.sgn);
        ans.red = dec_frac1(ans.num, ans.den) || dec_frac0(ans.num, ans.den) || dec_is_zero(ans.den);
        return ans;
    }
    net_decimal &operator/=(const net_decimal &divr) {
        *this = *this / divr;
        return *this;
    }
    callback_dec_arg friend arg operator/=(arg &divd, const net_decimal &divr) { return divd /= divr.to_float(); }

    friend net_decimal operator%(const net_decimal &divd, const net_decimal &divr) {
        net_decimal ans;
        dec_frac_rem(ans.num, ans.den, divd.num, divd.den, divr.num, divr.den);
        if (divr.modulus_mode && divr.sgn) ans += divr;
        if (divd.sgn) ans.sgn = !ans.sgn;
        return ans;
    }
    net_decimal &operator%=(const net_decimal &divr) {
        *this = *this % divr;
        return *this;
    }
    callback_dec_arg friend arg operator%=(arg &divd, const net_decimal &divr) {
        auto divr_int = divr.to_int();
        divd         %= divr_int;
        if (divr.modulus_mode && divr.sgn) divd += divr_int;
        if (divd < 0) divd *= -1;
        return divd;
    }

    friend net_decimal operator<<(const net_decimal &src, const net_decimal &bit) {
        auto ans = src;
        ans <<= bit;
        return ans;
    }
    net_decimal &operator<<=(const net_decimal &bit) {
        net_decimal_data bit_cnt;
        if (!dec_frac_bit_verify(bit_cnt, bit.num, bit.den)) return *this;
        reduct();
        if (dec_is_zero(den)) dec_bit_lsh(num, bit_cnt.it[0], binary_bit_set);
        return *this;
    }
    callback_dec_arg friend arg operator<<=(arg &src, const net_decimal &bit) {
        net_decimal_data bit_cnt;
        if (!dec_frac_bit_verify(bit_cnt, bit.num, bit.den)) return src;
        if (dec_is_zero(den)) src <<= bit_cnt.it[0];
        return src;
    }

    friend net_decimal operator>>(const net_decimal &src, const net_decimal &bit) {
        auto ans = src;
        ans >>= bit;
        return ans;
    }
    net_decimal &operator>>=(const net_decimal &bit) {
        net_decimal_data bit_cnt;
        if (!dec_frac_bit_verify(bit_cnt, bit.num, bit.den)) return *this;
        reduct();
        if (dec_is_zero(den)) dec_bit_rsh(num, bit_cnt.it[0], binary_bit_set);
        return *this;
    }
    callback_dec_arg friend arg operator>>=(arg &src, const net_decimal &bit) {
        net_decimal_data bit_cnt;
        if (!dec_frac_bit_verify(bit_cnt, bit.num, bit.den)) return src;
        if (dec_is_zero(den)) src >>= bit_cnt.it[0];
        return src;
    }

    friend net_decimal operator&(const net_decimal &fst, const net_decimal &snd) { return binary_operator<NEUNET_DEC_BIN_AND>(fst, snd); }
    net_decimal &operator&=(const net_decimal &snd) {
        *this = *this & snd;
        return *this;
    }
    callback_dec_arg friend arg operator&=(arg &fst, const net_decimal &snd) {
        if (dec_is_zero(snd.den)) {
            if (snd.num.ft.length) return fst;
            return fst &= snd.num.it[0];
        }
        auto tmp = snd.num,
             opt = dec_rem(tmp, snd.den);
        if (dec_is_zero(tmp) && !tmp.ft.length) return fst &= opt.it[0];
        return fst;
    }

    friend net_decimal operator|(const net_decimal &fst, const net_decimal &snd) { return binary_operator<NEUNET_DEC_BIN_OR>(fst, snd); }
    net_decimal &operator|=(const net_decimal &snd) {
        *this = *this | snd;
        return *this;
    }
    callback_dec_arg friend arg operator|=(arg &fst, const net_decimal &snd) {
        if (dec_is_zero(snd.den)) {
            if (snd.num.ft.length) return fst;
            return fst |= snd.num.it[0];
        }
        auto tmp = snd.num,
             opt = dec_rem(tmp, snd.den);
        if (dec_is_zero(tmp) && !tmp.ft.length) return fst |= opt.it[0];
        return fst;
    }

    friend net_decimal operator^(const net_decimal &fst, const net_decimal &snd) { return binary_operator<NEUNET_DEC_BIN_XOR>(fst, snd); }
    net_decimal &operator^=(const net_decimal &snd) {
        *this = *this ^ snd;
        return *this;
    }
    callback_dec_arg friend arg operator^=(arg &fst, const net_decimal &snd) {
        if (dec_is_zero(snd.den)) {
            if (snd.num.ft.length) return fst;
            return fst ^= snd.num.it[0];
        }
        auto tmp = snd.num,
             opt = dec_rem(tmp, snd.den);
        if (dec_is_zero(tmp) && !tmp.ft.length) return fst ^= opt.it[0];
        return fst;
    }

    net_decimal operator~() const {
        net_decimal ans;
        if (dec_is_zero(den)) {
            if (num.ft.length) return *this;
            ans.num = num;
            dec_bit_not(ans.num, false, binary_bit_set, sgn);
            return ans;
        }
        auto rm = num;
        ans.num = dec_rem(rm, den);
        if (dec_is_zero(rm)) dec_bit_not(ans.num, false, binary_bit_set, sgn);
        else ans = *this;
        return ans;
    }

    friend std::ostream &operator<<(std::ostream &os, const net_decimal &src) {
        if (src.sgn) os << '-';
        if (dec_is_zero(src.den)) return os << src.num;
        return os << dec_div(src.num, src.den, src.division_precision);
    }

protected:

public:
};

net_decimal operator""_d(const char *src, uint64_t len) { return net_decimal(std::string(src)); }
net_decimal operator""_d(long double src) { return net_decimal(src); }
net_decimal operator""_d(uint64_t src) { return net_decimal(src); }

NEUNET_END

_STD_BEGIN

neunet::net_decimal abs(const neunet::net_decimal &_Xx) { return _Xx.abs(); }

// neunet::net_decimal pow(const neunet::net_decimal &_Xx, const neunet::net_decimal &_Yx) { return neunet::net_decimal::dec_pow(_Xx, _Yx); }

// neunet::net_decimal exp(const neunet::net_decimal &_Xx) { return neunet::net_decimal::dec_exp(_Xx); }

// neunet::net_decimal log(const neunet::net_decimal &_Xx) { return neunet::net_decimal::dec_ln(_Xx); }

// neunet::net_decimal sin(const neunet::net_decimal &_Xx) { return neunet::net_decimal::dec_sin(_Xx); }

// neunet::net_decimal cos(const neunet::net_decimal &_Xx) { return neunet::net_decimal::dec_cos(_Xx); }

_STD_END