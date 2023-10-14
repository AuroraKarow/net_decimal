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

/** 
 * @brief Left shift
 * @param src   [IO]    source
 * @param bit   [In]    shift count
 * @param sta   [In]    bit stability
 */
void dec_bit_lsh(net_decimal_data &src, uint64_t bit = 1, bool sta = false){
    if (src.ft.length || !src.it.length) return;
    uint64_t b = bit % NEUNET_DEC_DIG_MAX,
             B = bit / NEUNET_DEC_DIG_MAX,
             l = src.it.length,
             T = 0,
             t;
    if (!sta && (B || src.it[l - 1] >= (NEUNET_DEC_SEG_MAX >> b))) src.it.init(l + (bit / 63) + 1);
    if (src.it[l - 1] >= (NEUNET_DEC_SEG_MAX >> b)) l++;
    // std::cout << "S" << std::endl;
    if (b){
        uint64_t bas = NEUNET_DEC_SEG_MAX >> b;
        for(uint64_t m = 0; m < l; m++) {
            t = src.it[m] / bas;
            src.it[m] %= bas;
            src.it[m] <<= b;
            src.it[m] += T;
            T = t;
        }
    }
    // std::cout << "A" << std::endl;
    uint64_t Bas = NEUNET_DEC_SEG_MAX >> NEUNET_DEC_DIG_MAX;
    while (B--){
        T = 0;
        for(uint64_t m = 0; m < src.it.length; m++) {
            t = src.it[m] / Bas;
            src.it[m] %= Bas;
            src.it[m] <<= NEUNET_DEC_DIG_MAX;
            src.it[m] += T;
            T = t;
        }
    }
    l = src.it.length;
    // std::cout << "B" << std::endl;
    if (!sta) {
        while(!src.it[--l]) continue;
        src.it = src.it.sub_set(0, l);
    }
}

/**
 * @brief Right shift
 * @param src   [IO]    source
 * @param bit   [In]    shift count
 * @param sta   [In]    bit stability
 */
void dec_bit_rsh(net_decimal_data &src, uint64_t bit = 1, bool sta = false){
    if (src.ft.length || !src.it.length) return;
    uint64_t b = bit % NEUNET_DEC_DIG_MAX,
             B = bit / NEUNET_DEC_DIG_MAX,
             l = src.it.length,
             T = 0,
             t;
    // std::cout << "S" << std::endl;
    if (b){
        uint64_t bas = NEUNET_DEC_SEG_MAX >> b,
                 div = 1 << b;
        for(uint64_t m = l; m > 0; m--) {
            t = (src.it[m - 1] % div) * bas;
            src.it[m - 1] >>= b;
            src.it[m - 1] += T;
            T = t;
        }
    }
    // std::cout << "A" << std::endl;
    uint64_t Bas = NEUNET_DEC_SEG_MAX >> NEUNET_DEC_DIG_MAX,
             Div = 1 << NEUNET_DEC_DIG_MAX;
    while (B--){
        T = 0;
        for(uint64_t m = l; m > 0; m--) {
            t = (src.it[m - 1] % Div) * Bas;
            src.it[m - 1] >>= NEUNET_DEC_DIG_MAX;
            src.it[m - 1] += T;
            T = t;
        }
    }
    // std::cout << "B" << std::endl;
    if (!sta) {
        while(!src.it[--l]) continue;
        src.it = src.it.sub_set(0, l);
    }
}

void dec_bit_lsh1(net_decimal_data &src, bool sta = false){
    if (src.ft.length || !src.it.length) return;
    if (!sta && src.it[src.it.length - 1] >= NEUNET_DEC_BIT_BAS) src.it.init(src.it.length + 1);
    int t = 0;
    for(uint64_t m = 0; m < src.it.length; m++) if(src.it[m] >= NEUNET_DEC_BIT_BAS) {
        src.it[m] -= NEUNET_DEC_BIT_BAS;
        src.it[m] <<= 1;
        src.it[m] += t;
        t = 1;
    } else {
        src.it[m] <<= 1;
        src.it[m] += t;
        t = 0;
    }
}

void dec_bit_rsh1(net_decimal_data &src, bool sta = false){
    if (src.ft.length || !src.it.length) return;
    //if (!(src.it.length - 1)) src.it[0] >>= 1; return;
    uint64_t t = 0;
    for(uint64_t m = src.it.length; m > 0; m--) if(src.it[m - 1] % 2 == 1){ 
        src.it[m - 1] >>= 1;
        src.it[m - 1] += t;
        t = NEUNET_DEC_BIT_BAS;
    } else{
        src.it[m - 1] >>= 1;
        src.it[m - 1] += t;
        t = 0;
    }
    if (!sta && !src.it[src.it.length - 1]) src.it.length - 1 ? src.it[src.it.length - 2] < NEUNET_DEC_BIT_TOP ? src.it.init(src.it.length - 1, true) : src.it.init(src.it.length, true) : src.it.init(src.it.length - 1, true);
}

// kp
net_decimal_data dec_bit_k0(bool sgn = false) { return dec_init(sgn, 0); }
// k, k1
net_decimal_data dec_bit_k1(bool sgn = false) { return dec_init(sgn, 1); }

uint64_t dec_bit_cnt(net_decimal_data &src){
    uint64_t b = 0;
    auto K = src,
         k = dec_bit_k0();
    while(dec_comp(k, K) != NEUNET_DEC_CMP_EQL){ 
        dec_bit_rsh1(K);
        b++;
    }
    return b;
}

/**
 * @brief Binary not
 * @param src       [IO]    source
 * @param comple    [In]    one's complement
 * @param bit       [In]    lower bit count, minimum is the bit count of source value
 * @param sgn       [In]    Sign of source
 */
void dec_bit_not(net_decimal_data &src, bool comple = false, uint64_t bit = 0, bool sgn = false) {
    auto k0 = dec_bit_k0(),
         k1 = dec_bit_k1(),
         k  = k0,
         K  = k1;
    auto b  = 0;
    auto k_sgn = false;
    if (dec_comp(src, k0) == NEUNET_DEC_CMP_EQL && !bit) {src = std::move(k1); return;}
    if (comple && !sgn) return;
    while(dec_comp(src, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (src.it.length || !(src.it.length - 1) && src.it[0] != 0) { if (!(src.it[0] % 2)) k = dec_add(k_sgn, k, false, k1, false); }
        else if (!src.it.length) k = dec_add(k_sgn, k, false, k1, false);
        dec_bit_lsh1(k1);
        dec_bit_rsh1(src);
        b++;
    }
    if (comple) {
        k = dec_add(k_sgn, k, false, K, false);
        if (sgn) { if (dec_comp(k, k1) != NEUNET_DEC_CMP_LES) { dec_bit_lsh1(k1); } k = dec_add(k_sgn, k, false, k1, false); }
    }
    else if (!sgn) k = dec_add(k_sgn, k, false, k1, false);
    /*
    while(dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
        if(!(D1.it[0] % 2)){ k = dec_add(sgn, k, sgn, kp, sgn); }
        dec_bit_lsh_one(kp);
        dec_bit_rsh_one(D1);
    }
    */
    src = std::move(k);
}

/**
 * @brief Binary not
 * @param fst   [In]    first source
 * @param snd   [In]    second source
 * @param sgn1  [In]    first source sign
 * @param sgn2  [In]    second source sign
 * @param bit   [In]    lower bit count, minimum is the bit count of sources value
 */
net_decimal_data dec_bit_and(const net_decimal_data &fst, const net_decimal_data &snd, bool sgn1 = false, bool sgn2 = false, uint64_t bit = 0) {
    auto D  = fst,
         d  = snd,
         k0 = dec_bit_k0(),
         k1 = dec_bit_k1(),
         k  = k0;
    auto b  = 0;
    auto k_sgn = false;
    if (!bit) {
        auto a = dec_bit_cnt(D),
             b = dec_bit_cnt(d);
        bit = std::max(a, b);
    }
    if (sgn1) dec_bit_not(D, true, bit, true);
    if (sgn2) dec_bit_not(d, true, bit, true);
    while(dec_comp(d, k0) != NEUNET_DEC_CMP_EQL || dec_comp(D, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (!(!D.it.length || !d.it.length)) if (D.it[0] % 2 && d.it[0] % 2) k = dec_add(k_sgn, k, false, k1, false);
        dec_bit_lsh1(k1);
        dec_bit_rsh1(D);
        dec_bit_rsh1(d);
        b++;
    }
    /*
    while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
        if(D1.it[0] % 2 && d1.it[0] % 2){ k = dec_add(sgn, k, sgn, kp, sgn); }
        dec_bit_lsh_one(kp);
        dec_bit_rsh_one(D1);
        dec_bit_rsh_one(d1);
    }
    */
    return k;
}

/**
 * @brief Binary not
 * @param fst   [In]    first source
 * @param snd   [In]    second source
 * @param sgn1  [In]    first source sign
 * @param sgn2  [In]    second source sign
 * @param bit   [In]    lower bit count, minimum is the bit count of sources value
 */
net_decimal_data dec_bit_or(const net_decimal_data &fst, const net_decimal_data &snd, bool sgn1 = false, bool sgn2 = false, uint64_t bit = 0) {
    auto D  = fst,
         d  = snd,
         k0 = dec_bit_k0(),
         k1 = dec_bit_k1(),
         k  = k0;
    auto b  = 0;
    auto k_sgn = false;
    if (!bit) {
        auto a = dec_bit_cnt(D),
             b = dec_bit_cnt(d);
        bit = a > b ? a : b;
    }
    if (sgn1) dec_bit_not(D, true, bit, true);
    if (sgn2) dec_bit_not(d, true, bit, true);
    while(dec_comp(d, k0) != NEUNET_DEC_CMP_EQL || dec_comp(D, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (!(!D.it.length || !d.it.length)) { if((D.it[0] % 2) || (d.it[0] % 2)) k = dec_add(k_sgn, k, true, k1, true); }
        else if (!D.it.length) { if (d.it[0] % 2) k = dec_add(k_sgn, k, false, k1, false); }
        else if (!d.it.length) { if (D.it[0] % 2) k = dec_add(k_sgn, k, false, k1, false); }
        dec_bit_lsh1(k1);
        dec_bit_rsh1(D);
        dec_bit_rsh1(d);
        b++;
    }
    /*
    while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
        if(!(!(D1.it[0] % 2) && !(d1.it[0] % 2))){ k = dec_add(sgn, k, sgn, kp, sgn); }
        dec_bit_lsh_one(kp);
        dec_bit_rsh_one(D1);
        dec_bit_rsh_one(d1);
    }
    */
    return k;
}

/**
 * @brief Binary not
 * @param fst   [In]    first source
 * @param snd   [In]    second source
 * @param sgn1  [In]    first source sign
 * @param sgn2  [In]    second source sign
 * @param bit   [In]    lower bit count, minimum is the bit count of sources value
 */
net_decimal_data dec_bit_xor(const net_decimal_data &fst, const net_decimal_data &snd, bool sgn1 = false, bool sgn2 = false, uint64_t bit = 0) {
    auto k0 = dec_bit_k0(),
         k1 = dec_bit_k1(),
         k  = k0,
         D  = fst,
         d  = snd;
    auto b  = 0;
    auto k_sgn = false;
    if (!bit) {
        auto a = dec_bit_cnt(D),
             b = dec_bit_cnt(d);
        bit = a > b ? a : b;
    }
    if (sgn1) dec_bit_not(D, true, bit, true);
    if (sgn2) dec_bit_not(d, true, bit, true);
    while(dec_comp(d, k0) != NEUNET_DEC_CMP_EQL || dec_comp(D, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (!(!D.it.length || !d.it.length)) { if((D.it[0] % 2 && !(d.it[0] % 2)) || (!(D.it[0] % 2) && d.it[0] % 2)) k = dec_add(k_sgn, k, true, k1, true); }
        else if (!D.it.length) { if (d.it[0] % 2) k = dec_add(k_sgn, k, false, k1, false); }
        else if (!d.it.length) { if (D.it[0] % 2) k = dec_add(k_sgn, k, false, k1, false); }
        dec_bit_lsh1(k1);
        dec_bit_rsh1(D);
        dec_bit_rsh1(d);
        b++;
    }
    /*
    while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
        if((D1.it[0] % 2 && !(d1.it[0] % 2)) || (!(D1.it[0] % 2) && d1.it[0] % 2)){ k = dec_add(sgn, k, sgn, kp, sgn); }
        dec_bit_lsh_one(kp);
        dec_bit_rsh_one(D1);
        dec_bit_rsh_one(d1);
    }
    */
    return k;
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
            d = dec_e10(s, precision + 2);
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
        d = dec_e10(x, prec + 2);
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
        d = dec_e10(u, prec + 2);
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
        d = dec_e10(u, prec + 2);
    } while (dec_comp(d, v) == NEUNET_DEC_CMP_GTR);
    return dec_div(ans_sgn, p, p_sgn, q, false, prec);
}



NEUNET_END