NEUNET_BEGIN

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

// NOT thread safety
class {
private:
    net_set<net_decimal_data> prec;
    
public:
    net_decimal_data get(uint64_t precision) {
        if (!precision) return dec_e10(0);
        --precision;
        if (precision >= prec.length) {
            auto len = prec.length;
            if (!len) len = 128;
            while (len < precision) len <<= 1;
            prec.init(len);
        }
        if (dec_is_zero(prec[precision])) prec[precision] = dec_e10(precision + 1);
        return prec[precision];
    }

} net_decimal_prec;

// NOT thread safety
class {
private:
    net_decimal_data b,
                     c,
                     s,
                     o,
                     p,
                     q,
                     a;

    bool need_init = true;

    uint64_t prec = 0;

public:
    net_decimal_data get(uint64_t precision) {
        if (precision <= prec) return a;
        prec = precision;
        if (need_init) {
            need_init = false;
            auto sgn  = false;
            b = dec_init(sgn, 1);
            c = dec_init(sgn, 2);
            s = dec_init(sgn, "0.6");
            o = dec_init(sgn, "0.36");
            q = b;
        }
        net_decimal_data d;
        do {
            if (dec_comp(b, q) == NEUNET_DEC_CMP_EQL) p = dec_add(p, s);
            else {
                p = dec_add(dec_mul(b, p), dec_mul(s, q));
                q = dec_mul(q, b);
            }
            s = dec_mul(s, o);
            b = dec_add(b, c);
            d = dec_e10_mul(s, precision + 2);
        } while (dec_comp(d, b) == NEUNET_DEC_CMP_GTR);
        a = dec_mul(c, dec_div(p, q, prec));
        return a;
    }

} net_decimal_ln4;

// NOT thread safety
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

    void run(uint64_t precision) { while (prec < precision) {
        init();
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
    } }

public:
    net_decimal_data get(uint64_t precision) {
        run(precision);
        return dec_div(p, q, precision);
    }

    const net_decimal_data &num(uint64_t precision) {
        run(precision);
        return p;
    }

    net_decimal_data num_2pi(uint64_t precision) {
        run(precision);
        return dec_mul(dec_2, p);
    }

    const net_decimal_data &den(uint64_t precision) {
        run(precision);
        return q;        
    }

} net_decimal_pi;

// 0 < src <= 2
net_decimal_data dec_ln_2(bool &ans_sgn, const net_decimal_data &src, uint64_t prec) {
    auto p_sgn = false,
         x_sgn = false,
         c_sgn = false;
    auto b = dec_init(p_sgn, 1),
         p = dec_init(p_sgn, 0),
         i = b,
         q = b,
         c = b,
         x = dec_sub(x_sgn, src, false, b, false),
         o = x,
         d = p;
    auto o_sgn = x_sgn;
    do {
        if (dec_comp(b, q) == NEUNET_DEC_CMP_EQL) p = dec_add(p_sgn, p, p_sgn, x, x_sgn);
        else {
            auto fst_sgn = false,
                 snd_sgn = false;
            auto fst = dec_mul(fst_sgn, b, false, p, p_sgn),
                 snd = dec_mul(snd_sgn, c, c_sgn, x, x_sgn);
            snd = dec_mul(snd_sgn, snd, snd_sgn, q, false);
            p   = dec_add(p_sgn, fst, fst_sgn, snd, snd_sgn);
            q   = dec_mul(q, b);
        }
        x = dec_mul(x_sgn, x, x_sgn, o, o_sgn);
        b = dec_add(b, i);
        c_sgn = !c_sgn;
        d = dec_e10_mul(x, prec + 2);
    } while (dec_comp(d, b) == NEUNET_DEC_CMP_GTR);
    return dec_div(ans_sgn, p, p_sgn, q, false, prec);
}

net_decimal_data dec_ln(bool &ans_sgn, const net_decimal_data &src, uint64_t prec) {
    net_decimal_data cnt,
                     one = dec_init(ans_sgn, 1);
    auto cmp = dec_comp(src, one);
    if (cmp == NEUNET_DEC_CMP_EQL) return cnt;
    if (cmp == NEUNET_DEC_CMP_LES) return dec_ln_2(ans_sgn, src, prec);
    net_decimal_data base = src,
                     ofhs = dec_init(ans_sgn, 0.25);
    while (dec_comp(base, one) == NEUNET_DEC_CMP_GTR) {
        base = dec_mul(base, ofhs);
        cnt  = dec_add(cnt, one);
    }
    auto cnt_sgn = false;
    cnt  = dec_mul(cnt_sgn, cnt, false, net_decimal_ln4.get(prec), false);
    base = dec_ln_2(ans_sgn, base, prec);
    return dec_add(ans_sgn, base, ans_sgn, cnt, cnt_sgn);
}

// TODO waiting for optimizing ... 
net_decimal_data dec_exp(bool &ans_sgn, const net_decimal_data &src_num, const net_decimal_data &src_den,  bool src_sgn, uint64_t prec) {
    auto c = dec_init(ans_sgn, 1),
         p = c,
         q = c,
         o = c,
         u = src_num,
         v = src_den,
         d = dec_init(ans_sgn, 0);
    auto p_sgn = false,
         u_sgn = src_sgn,
         s_dec = dec_comp(src_den, d) == NEUNET_DEC_CMP_EQL;
    if (s_dec) v = c;
    do {
        if (dec_comp(q, v) == NEUNET_DEC_CMP_EQL) p = dec_add(p_sgn, p, p_sgn, u, u_sgn);
        else {
            auto fst_sgn = false,
                 snd_sgn = false;
            auto fst = dec_mul(fst_sgn, v, false, p, p_sgn),
                 snd = dec_mul(snd_sgn, u, u_sgn, q, false);
            p = dec_add(p_sgn, fst, fst_sgn, snd, snd_sgn);
            q = dec_mul(q, v);
        }
        c = dec_add(c, o);
        u = dec_mul(u_sgn, u, u_sgn, src_num, src_sgn);
        if (s_dec) v = dec_mul(v, c);
        else v = dec_mul(v, dec_mul(c, src_den));
        d = dec_e10_mul(u, prec + 2);
    } while (dec_comp(d, v) == NEUNET_DEC_CMP_GTR);
    return dec_div(ans_sgn, p, p_sgn, q, false, prec);
}

net_decimal_data dec_sin_cos(bool &ans_sgn, const net_decimal_data &src_num, const net_decimal_data &src_den, bool src_sgn, net_decimal_data &u, net_decimal_data &v, net_decimal_data &b, uint64_t prec) {
    auto p_sgn = false,
         u_sgn = false,
         m_sgn = false,
         c_sgn = false;
    auto c = dec_init(ans_sgn, 1),
         q = c,
         p = dec_init(ans_sgn, 0),
         m = dec_mul(m_sgn, src_num, src_sgn, src_num, src_sgn),
         n = c,
         d = p;
    auto frac = dec_comp(src_den, d) != NEUNET_DEC_CMP_EQL;
    if (frac) n = dec_mul(src_den, src_den);
    else v = c;
    do {
        auto t_sgn = false;
        auto t     = dec_mul(t_sgn, c, c_sgn, u, u_sgn);
        if (dec_comp(q, v) == NEUNET_DEC_CMP_EQL) p = dec_add(p_sgn, p, p_sgn, t, t_sgn);
        else {
            p = dec_mul(p_sgn, v, false, p, p_sgn);
            t = dec_mul(t_sgn, t, t_sgn, q, false);
            p = dec_add(p_sgn, p, p_sgn, t, t_sgn);
            q = dec_mul(q, v);
        }
        u = dec_mul(u_sgn, u, u_sgn, m, m_sgn);
        if (frac) v = dec_mul(v, n);
        for (auto i = 0; i < 2; ++i) {
            b = dec_add(b, c);
            v = dec_mul(v, b);
        }
        c_sgn = !c_sgn;
        d = dec_e10_mul(u, prec + 2);
    } while (dec_comp(d, v) == NEUNET_DEC_CMP_GTR);
    return dec_div(ans_sgn, p, p_sgn, q, false, prec);
}

NEUNET_END