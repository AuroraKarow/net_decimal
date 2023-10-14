NEUNET_BEGIN

class {
private:
    bool on  = false,
         sgn = false;

    std::atomic_uint64_t prec = 0;

    net_decimal_data p, q,
                     // +2
                     a,
                     // +4
                     b, 
                     // +8
                     c, d,
                     o, s,
                     dec_2,
                     dec_4,
                     dec_8;
    
    void init() {
        if (on) return;
        p = dec_init(sgn, 0);
        q = dec_init(sgn, 1);
        a = dec_init(sgn, 0.25);
        b = dec_init(sgn, 2);
        c = dec_init(sgn, 5);
        d = dec_init(sgn, 6);
        o = dec_init(sgn, 0.0625);
        s = q;
        dec_2 = dec_init(sgn, 2);
        dec_4 = dec_init(sgn, 4);
        dec_8 = dec_init(sgn, 8);
        on = true;
    }

public:
    net_decimal_data get(uint64_t precision) {
        init();
        while (prec < precision) {
            
            auto u = dec_init(sgn, 1),
                 v = a;
            u = dec_add(dec_mul(b, u), v, true);
            v = dec_mul(v, b);
            u = dec_add(dec_mul(c, u), v, true);
            v = dec_mul(v, c);
            u = dec_add(dec_mul(d, u), v, true);
            v = dec_mul(v, d);
            u = dec_mul(u, s);
            if (dec_comp(q, v) == NEUNET_DEC_CMP_EQL) p = dec_add(p, u);
            else {
                p = dec_add(dec_mul(v, p), dec_mul(u, q));
                q = dec_mul(v, q);
            }
            s = dec_mul(s, o);
            a = dec_add(a, dec_2);
            b = dec_add(b, dec_4);
            c = dec_add(c, dec_8);
            d = dec_add(d, dec_8);
            ++prec;
        }
        return dec_div(p, q, precision);
    }

} net_decimal_pi;

class net_decimal {
public:
    bool modulus_mode    = false,
         fraction_reduct = false;

    uint64_t division_precision = 32;

protected:
    void value_assign(const net_decimal &src) {
        sgn                = src.sgn;
        dec                = src.dec;
        modulus_mode       = src.modulus_mode;
        fraction_reduct    = src.fraction_reduct;
        division_precision = src.division_precision;
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

    bool is_integer() {
        if (dec) return !base.ft.length;
        dec_reduct(frac);
        if (!dec_is_one(frac.den)) return false;
        dec  = true;
        base = frac.num;
        return true;
    }

    void frac_init(net_decimal &u, net_decimal &v) {
        if (dec) {
            u = *this;
            v = 1;
        } else {
            u.sgn  = sgn;
            u.base = frac.num;
            v.base = frac.den;
        }
    }

    static net_decimal ln_4(uint64_t prec, bool show_iter = false) {
        net_decimal b {1},
                    c {2},
                    s {"0.6"},
                    o {"0.36"},
                    p {0},
                    q {b},
                    m {p},
                    n {p};
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
            n.base = dec_div(n.sgn, p.base, p.sgn, q.base, q.sgn, prec);
        } while (m != n);
        return 2 * n;
    }
    // ln(x + 1) = (-1 < x <= 1)
    net_decimal ln_2(uint64_t prec, bool show_iter = false) {
        net_decimal b {1},
                    p {0},
                    q {b},
                    x {*this - 1},
                    o {x},
                    c {b},
                    m {p},
                    n {p};
        do {
            if (show_iter) std::cout << n << std::endl;
            m = std::move(n);
            if (b == q) p += x;
            else {
                p  = b * p + c * x * q;
                q *= b;
            }
            x *= o;
            ++b;
            n.base = dec_div(n.sgn, p.base, p.sgn, q.base, q.sgn, prec);
            c.sgn  = !c.sgn;
        } while (m != n);
        return n;
    }
    net_decimal ln_proto(uint64_t prec, bool show_iter = false) {
        net_decimal b {1},
                    p {0},
                    q {b},
                    o {p},
                    u {*this - 1},
                    v {*this + 1},
                    h {u * u},
                    s {v * v},
                    m {p},
                    n {p};
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
            n.base = dec_div(n.sgn, p.base, p.sgn, q.base, q.sgn, prec);
        } while (m != n);
        return 2 * n;
    }

    net_decimal sin_cos_proto(net_decimal &u, net_decimal &v, net_decimal &b) {
        net_decimal c {1},
                    q {c},
                    p {0},
                    g {p},
                    h {p},
                    m {p},
                    n {p},
                    o {*this * (*this)};
        if (!dec) {
            g.base = o.frac.num;
            h.base = o.frac.den;
        }
        do {
            m = std::move(n);
            if (q == v) p += c * u;
            else {
                p  = v * p + c * q * u;
                q *= v;
            }
            if (dec) u *= o;
            else {
                u *= g;
                v *= h;
            }
            n.base = dec_div(n.sgn, p.base, p.sgn, q.base, q.sgn, division_precision);
            for (auto i = 0; i < 2; ++i) v *= (++b);
            c.sgn = !c.sgn;
        } while (m != n);
        return n;
    }

public:
    net_decimal() {}
    net_decimal(const net_decimal &src) { value_copy(src); }
    net_decimal(net_decimal &&src) { value_move(std::move(src)); }
    callback_arg net_decimal(const arg &src) {
        static_assert(neunet_dec_init_expr);
        base = dec_init(sgn, src);
    }

    std::string to_string() {
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
        division_precision = 32;
        fraction_reduct    = false;
        modulus_mode       = false;
        dec                = true;
        sgn                = false;
    }

    ~net_decimal() { reset(); }

    net_decimal abs() const {
        auto tmp (*this);
        tmp.sgn = false;
        return tmp;
    }

    net_decimal ln() {
        if (is_zero() || is_one() || sgn) return {};
        if (!dec) {
            net_decimal num, den;
            num.base = frac.num;
            den.base = frac.den;
            return num.ln() / den.ln();
        }
        net_decimal cnt  {0},
                    base {*this};
        while (base > 1) {
            base *= 0.25;
            ++cnt;
        }
        return base.ln_2(division_precision) + cnt * ln_4(division_precision);
    }

    net_decimal exp() {
        net_decimal c {1},
                    p {c},
                    q {c},
                    v {0},
                    u {v},
                    a {v},
                    b {v},
                    m {p},
                    n {p};
        frac_init(u, v);
        frac_init(a, b);
        do {
            m = std::move(n);
            if (q == v) p += u;
            else {
                p  = v * p + u * q;
                q *= v;
            }
            ++c;
            if (dec) {
                u *= *this;
                v *= c;
            } else {
                u *= a;
                v *= c * b;
            }
            n.base = dec_div(n.sgn, p.base, p.sgn, q.base, q.sgn, division_precision);
        } while (m != n);
        return n;
    }

    net_decimal sin() {
        net_decimal  b {1},
                     u {0},
                     v {u};
        frac_init(u, v);
        return sin_cos_proto(u, v, b);
    }

    net_decimal cos() {
        net_decimal b {0},
                    u {1},
                    v {u};
        return sin_cos_proto(u, v, b);
    }

    callback_dec net_decimal pow(dec_t &&times) {
        neunet_type_if (neunet_dec_num(dec_t)) return pow(net_decimal{times});
        neunet_type_else
        if (is_zero()) return 0.;
        if ((!sgn && is_one()) || times.is_zero()) return 1;
        if (times.is_one()) {
            if (times.sgn) return 1 / (*this);
            return *this;
        }
        if (!sgn) return (times * ln()).exp();
        net_decimal piv {0};
        piv.base = net_decimal_pi.get(division_precision + 1);
        net_decimal ans {0},
                    pia {piv * times},
                    prd {ans};
        net_decimal res[NEUNET_CACHE_LEN];
        uint64_t    res_cnt {0},
                    cnt {1};
        do {
            prd = cnt * pia;
            ans = prd.sin();
            // dec_truncate(ans.base, division_precision - 1);
            // std::cout << ans << '\n' << std::endl;
            if (ans.is_zero()) {
                prd = prd.cos();
                dec_truncate(prd.base, division_precision);
                // std::cout << prd << std::endl;
                if (sgn) return abs().pow(times) * prd;
                return pow(times) * prd;
            }
            res[res_cnt++] = (ans);
            cnt += 2;
        } while (res_cnt == 1 || ans != res[0]);
        return ans;
        neunet_type_endif
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
        if (!dec_is_valid(fst.frac)) fst.frac = dec_init(fst.base, fst.fraction_reduct);
        if (!dec_is_valid(snd.frac)) snd.frac = dec_init(snd.base, snd.fraction_reduct);
        ans.dec  = false;
        ans.frac = dec_add(ans.sgn, fst.frac, fst.sgn, snd.frac, snd.sgn, false);
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
        if (!dec_is_valid(minu.frac)) minu.frac = dec_init(minu.base, minu.fraction_reduct);
        if (!dec_is_valid(subt.frac)) subt.frac = dec_init(subt.base, subt.fraction_reduct);
        ans.dec  = false;
        ans.frac = dec_sub(ans.sgn, minu.frac, minu.sgn, subt.frac, subt.sgn, false);
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
        if (fst.dec && snd.dec) {
            ans.base = dec_mul(ans.sgn, fst.base, fst.sgn, snd.base, snd.sgn);
            return ans;
        }
        if (!dec_is_valid(fst.frac)) fst.frac = dec_init(fst.base, fst.fraction_reduct);
        if (!dec_is_valid(snd.frac)) snd.frac = dec_init(snd.base, snd.fraction_reduct);
        ans.dec  = false;
        ans.frac = dec_mul(ans.sgn, fst.frac, fst.sgn, snd.frac, snd.sgn, false);
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
            ans.frac = dec_init(divd.base, divr.base, false);
            return ans;
        }
        if (!dec_is_valid(divd.frac)) divd.frac = dec_init(divd.base, divd.fraction_reduct);
        if (!dec_is_valid(divr.frac)) divr.frac = dec_init(divr.base, divr.fraction_reduct);
        ans.frac = dec_div(ans.sgn, divd.frac, divd.sgn, divr.frac, divr.sgn, false);
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
            if (divr.modulus_mode && divd.sgn != divr.sgn) ans.base = dec_add(ans.sgn, ans.base, false, divr.base, divr.sgn);
        } else {
            net_decimal_data divr_base;
            auto tmp = dec_float_part(divr_base, divr.frac);
            net_assert(!dec_is_zero(tmp), "net_decimal", "%", "Divisor could not be 0.");
            dec_rem(ans.base, divr_base);
            if (divr.modulus_mode && divd.sgn != divr.sgn) ans.base = dec_add(ans.sgn, ans.base, false, divr_base, divr.sgn);
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
        if (divr.modulus_mode && divd_sgn != divr.sgn) divd += divr_tmp;
    }

    callback_dec_s friend net_decimal operator<<(fst_dec_t &&src, snd_dec_t &&bit) {
        neunet_type_if (neunet_dec_num(fst_dec_t)) return net_decimal {src} << bit;
        neunet_type_elif (neunet_dec_num(snd_dec_t)) return src << net_decimal {bit};
        neunet_type_else
        net_assert(src.is_integer() && bit.is_integer(), "net_decimal", "<<", "Variables should be ineteger.");
        net_decimal ans = src;
        // TODO lsh operation
        // for (net_decimal i = 0; i < bit; ++i) dec_bit_lsh_one(ans.base);
        return ans;
        neunet_type_endif
    }
    callback_dec void operator<<=(dec_t &&src) {
        neunet_type_if (neunet_dec_num(dec_t)) *this = *this << net_decimal {src};
        neunet_type_else *this = *this << src;
        neunet_type_endif
    }
    callback_dec_arg_ref friend void operator<<=(arg &fst, net_decimal &snd) { fst <<= snd.to_integer(); }

    callback_dec_s friend net_decimal operator>>(fst_dec_t &&src, snd_dec_t &&bit) {
        neunet_type_if (neunet_dec_num(fst_dec_t)) return net_decimal {src} << bit;
        neunet_type_elif (neunet_dec_num(snd_dec_t)) return src << net_decimal {bit};
        neunet_type_else
        net_assert(src.is_integer() && bit.is_integer(), "net_decimal", ">>", "Variables should be ineteger.");
        net_decimal ans = src;
        // TODO rsh operation
        // for (net_decimal i = 0; i < bit; ++i) dec_bit_rsh_one(ans.base);
        return ans;
        neunet_type_endif
    }
    callback_dec void operator>>=(dec_t &&src) {
        neunet_type_if (neunet_dec_num(dec_t)) *this = *this >> net_decimal {src};
        neunet_type_else *this = *this >> src;
        neunet_type_endif
    }
    callback_dec_arg_ref friend void operator>>=(arg &fst, net_decimal &snd) { fst >>= snd.to_integer(); }

    friend std::ostream &operator<<(std::ostream &os, const net_decimal &src) {
        if (src.sgn) os << '-';
        if (src.dec) {
            os << src.base;
            return os;
        }
        auto tmp = dec_div(src.frac.num, src.frac.den, src.division_precision);
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

neunet::net_decimal pow(neunet::net_decimal &_Xx, neunet::net_decimal &_Yx) { return _Xx.pow(_Yx); }

neunet::net_decimal exp(neunet::net_decimal &_Xx) { return _Xx.exp(); }

neunet::net_decimal log(neunet::net_decimal &_Xx) { return _Xx.ln(); }

neunet::net_decimal sin(neunet::net_decimal &_Xx) { return _Xx.sin(); }

neunet::net_decimal cos(neunet::net_decimal &_Xx) { return _Xx.cos(); }

_STD_END
