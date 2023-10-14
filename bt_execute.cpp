#include "iostream"
#include "net_decimal"
#include "bitset"
#include "math.h"
#include "chrono"

using namespace neunet;
#define NEUNET_CHRONO_TIME_POINT    std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()).time_since_epoch().count()

void dec_bit_lt_one(net_decimal_data &D, bool sta = false){
    if (D.ft.length || !D.it.length) return;
    if (!sta && D.it[D.it.length - 1] >= NEUNET_DEC_BIT_BAS) D.it.init(D.it.length + 1);
    int t = 0;
    for(uint64_t m = 0; m < D.it.length; m++) if(D.it[m] >= NEUNET_DEC_BIT_BAS) {
        D.it[m] -= NEUNET_DEC_BIT_BAS;
        D.it[m] <<= 1;
        D.it[m] += t;
        t = 1;
    } else {
        D.it[m] <<= 1;
        D.it[m] += t;
        t = 0;
    }
}

void dec_bit_rt_one(net_decimal_data &D, bool sta = false){
    if (D.ft.length || !D.it.length) return;
    if (!(D.it.length - 1)) D.it[0] >>= 1; return;
    uint64_t t = 0;
    for(uint64_t m = D.it.length; m > 0; m--) if(D.it[m - 1] % 2 == 1){ 
        D.it[m - 1] >>= 1;
        D.it[m - 1] += t;
        t = NEUNET_DEC_BIT_BAS;
    } else{
        D.it[m - 1] >>= 1;
        D.it[m - 1] += t;
        t = 0;
    }
    if (!sta && !D.it[D.it.length - 1]) D.it.length - 1 ? D.it[D.it.length - 2] < NEUNET_DEC_BIT_TOP ? D.it.init(D.it.length - 1, true) : D.it.init(D.it.length, true) : D.it.init(D.it.length - 1, true);
}

// kp
net_decimal_data dec_bit_k0(bool sgn = false) { return dec_init(sgn, 0); }
// k, k1
net_decimal_data dec_bit_k1(bool sgn = false) { return dec_init(sgn, 1); }

uint64_t dec_bit_cnt(net_decimal_data &D){
    uint64_t b = 0;
    auto K = D,
         k = dec_bit_k0();
    while(dec_comp(k, K) != NEUNET_DEC_CMP_EQL){ 
        dec_bit_rt_one(K);
        b++;
    }
    return b;
}

void dec_bit_not(net_decimal_data &src, bool comple = false, uint64_t bit = 0, bool sgn = false) {
    auto k0 = dec_bit_k0(),
         k1 = dec_bit_k1(),
         k  = k0,
         K  = k1;
    auto b  = 0;
    auto k_sgn = false;
    if (dec_comp(src, k0) == NEUNET_DEC_CMP_EQL && !bit) {src = std::move(k1); return;}
    if (comple && sgn) return;
    while(dec_comp(src, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if(dec_dig_cnt(src, true) || !(src.it.length - 1) && src.it[0] != 0) { if (!(src.it[0] % 2)) k = dec_add(k_sgn, k, true, k1, true); }
        else if (!dec_dig_cnt(src, true)) k = dec_add(k_sgn, k, true, k1, true);
        dec_bit_lt_one(k1);
        dec_bit_rt_one(src);
        b++;
    }
    if (comple) {
        k = dec_add(k_sgn, k, true, K, true);
        if (!sgn) { if (dec_comp(k, k1) != NEUNET_DEC_CMP_LES) { dec_bit_lt_one(k1); } k = dec_add(k_sgn, k, true, k1, true); }
    }
    else if (sgn) k = dec_add(k_sgn, k, true, k1, true);
    /*
    while(dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
        if(!(D1.it[0] % 2)){ k = dec_add(sgn, k, sgn, kp, sgn); }
        dec_bit_lsh_one(kp);
        dec_bit_rsh_one(D1);
    }
    */
    src = std::move(k);
}

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
        bit = a > b ? a : b;
    }
    if (!sgn1) dec_bit_not(D, true, bit);
    if (!sgn2) dec_bit_not(d, true, bit);
    while(dec_comp(d, k0) != NEUNET_DEC_CMP_EQL || dec_comp(D, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (!(!dec_dig_cnt(D, true) || !dec_dig_cnt(d, true))) if (D.it[0] % 2 && d.it[0] % 2) k = dec_add(k_sgn, k, true, k1, true);
        dec_bit_lt_one(k1);
        dec_bit_rt_one(D);
        dec_bit_rt_one(d);
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
    if (!sgn1) dec_bit_not(D, true, bit);
    if (!sgn2) dec_bit_not(d, true, bit);
    while(dec_comp(d, k0) != NEUNET_DEC_CMP_EQL || dec_comp(D, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (!(!dec_dig_cnt(D, true) || !dec_dig_cnt(d, true))) { if((D.it[0] % 2) || (d.it[0] % 2)) k = dec_add(k_sgn, k, true, k1, true); }
        else if (!dec_dig_cnt(D, true)) { if (d.it[0] % 2) k = dec_add(k_sgn, k, true, k1, true); }
        else if (!dec_dig_cnt(d, true)) { if (D.it[0] % 2) k = dec_add(k_sgn, k, true, k1, true); }
        dec_bit_lt_one(k1);
        dec_bit_rt_one(D);
        dec_bit_rt_one(d);
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

// net_decimal_data dec_rad_min(net_decimal_data &src, uint64_t bit = 0, bool sgn = true){
//     if (dec_comp(dec_bit_k0(), src) == NEUNET_DEC_CMP_EQL) return src;
//     auto S = src;
//     if (S.it[S.it.length - 1] >= NEUNET_DEC_BIT_TOP) S.it.init(S.it.length + 1);
//     if (!sgn) dec_bit_not(S, sgn, bit);
//     return S;
// }

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
    if (!sgn1) dec_bit_not(D, true, bit);
    if (!sgn2) dec_bit_not(d, true, bit);
    while(dec_comp(d, k0) != NEUNET_DEC_CMP_EQL || dec_comp(D, k0) != NEUNET_DEC_CMP_EQL || b < bit){
        if (!(!dec_dig_cnt(D, true) || !dec_dig_cnt(d, true))) { if((D.it[0] % 2 && !(d.it[0] % 2)) || (!(D.it[0] % 2) && d.it[0] % 2)) k = dec_add(k_sgn, k, true, k1, true); }
        else if (!dec_dig_cnt(D, true)) { if (d.it[0] % 2) k = dec_add(k_sgn, k, true, k1, true); }
        else if (!dec_dig_cnt(d, true)) { if (D.it[0] % 2) k = dec_add(k_sgn, k, true, k1, true); }
        dec_bit_lt_one(k1);
        dec_bit_rt_one(D);
        dec_bit_rt_one(d);
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

net_decimal_data dec_exec_bit(const int type, const net_decimal_data &D, const net_decimal_data &d){
    /*
    type:
        0 -> and
        1 -> or 
        2 -> not
        3 -> or not
    */
    bool sgn = false;
    auto k = dec_init(sgn, 0), k0 = dec_init(sgn, 0), kp = dec_init(sgn, 1), k1 = dec_init(sgn, 1), k2 = dec_init(sgn, 2), D1 = D, d1 = d;
    if(type == 0){
        while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
            if(D1.it[0] % 2 && d1.it[0] % 2){ k = dec_add(sgn, k, sgn, kp, sgn); }
            dec_bit_lt_one(kp);
            dec_bit_rt_one(D1);
            dec_bit_rt_one(d1);
        }
    }
    if(type == 1){
        while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
            if(!(!(D1.it[0] % 2) && !(d1.it[0] % 2))){ k = dec_add(sgn, k, sgn, kp, sgn); }
            dec_bit_lt_one(kp);
            dec_bit_rt_one(D1);
            dec_bit_rt_one(d1);
        }
    }
    if(type == 2){
        while(dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
            if(!(D1.it[0] % 2)){ k = dec_add(sgn, k, sgn, kp, sgn); }
            dec_bit_lt_one(kp);
            dec_bit_rt_one(D1);
        }
    }
    if(type == 3){
        while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
            if((D1.it[0] % 2 && !(d1.it[0] % 2)) || (!(D1.it[0] % 2) && d1.it[0] % 2)){ k = dec_add(sgn, k, sgn, kp, sgn); }
            dec_bit_lt_one(kp);
            dec_bit_rt_one(D1);
            dec_bit_rt_one(d1);
        }
    }
    return k;
}


net_decimal_data dec_bit_lt_one0(net_decimal_data &D, bool sgn = true){
    if (D.ft.length || !D.it.length) return D;
    auto d = D;
    dec_bit_not(d, true);
    dec_bit_lt_one(d);
    //if () D.it[D.it.length - 1] += k;
    return dec_bit_val(d);
}

net_decimal_data dec_bit_rt_one0(net_decimal_data &D, bool sgn = true){
    if (D.ft.length || !D.it.length) return D;
    auto d = D;
    dec_bit_not(d, true);
    dec_bit_rt_one(d);
    //if () D.it[D.it.length - 1] += k;
    return dec_bit_val(d);
}

uint64_t rev0(uint64_t &d){
    uint64_t dn = 0;
    for(uint64_t k = 0; k < 64; k++){ 
        uint64_t a = 1, b = 1;
        if(d & (a << k)){ dn += b << (63 - k); } 
    }
    d = dn;
    return d;
}

void dec_lt_move(net_decimal_data &D, int64_t B){ while(B--){ dec_bit_lt_one(D); } }

void dec_rt_move(net_decimal_data &D, int64_t B){ while(B--){ dec_bit_rt_one(D); } }



// void rev1(uint64_t &d){

// }

// void dec_rev_exec0(net_decimal_data &D, uint64_t *res_ds){
//     net_decimal_data d = D;
//     for(uint64_t i = 0; i < D.it.length; i++){ D.it[D.it.length - i - 1] = *(res_ds + d.it[i]); }
// }

// void dec_rev_exec1(net_decimal_data &D){
//     net_decimal_data d = D;
//     for(uint64_t i = 0; i < D.it.length; i++){ D.it[D.it.length - i - 1] = rev0(d.it[i]); }
// }
 
// void dec_rev_exec2(net_decimal_data &D){
//     for(uint64_t i = 0; i < D.it)
// }

net_decimal_data dec_bit_val(net_decimal_data &D){
    auto d = D,
         b = dec_bit_k1();
    bool sgn;
    auto bit = dec_bit_cnt(D);
    dec_lt_move(b, bit - 1);
    d = dec_add(sgn, d, false, dec_add(sgn, b, false, dec_bit_k1(), false), true);
    dec_bit_not(d, false, bit - 1);
    return d;
}


void main(){
    bool k_sgn = false;
    uint64_t k01 = 1;
    k01 <<= 63;
    uint64_t k02 = 1;
    k02 <<= 62;
    uint64_t k03 = 1;
    k03 <<= 55;
    uint64_t k = k01 + k02 + k03;
	std::cout << "k:" << std::endl;
	std::cout << k01 << std::endl;
	std::cout << k02 << std::endl;
	std::cout << k03 << std::endl;
	std::cout << k << std::endl;
    // auto k0 = dec_init(k_sgn, "190000000000000000019000000000000000001900000000000000000120000000000000000019");
    auto k0 = dec_init(k_sgn, "1"),
        k1 = dec_init(k_sgn, "10000000000000000000000000000000000000"),
        k2 = dec_init(k_sgn, "160"),
        k3 = dec_init(k_sgn, "206"),//11001110 / 306 / 100110010
        k4 = dec_init(k_sgn, "61");//00111101 / 451 / 111000011

    // auto start = NEUNET_CHRONO_TIME_POINT;
    // auto ds_rev = rev(64);
	// for(int a = 0;a < 99999;a++){ dec_rev_exec0(k0, ds_rev); }
	// auto end = NEUNET_CHRONO_TIME_POINT;
	// std::cout << "Time0: " << end - start << std::endl;
	// std::cout << dec_to_string(k_sgn, k0) << std::endl;
    auto K = k1, K1 = dec_bit_k1(), K2 = dec_bit_k1();
    dec_bit_lt_one(K1);
	std::cout << "<<: " << dec_to_string(k_sgn, K1) << std::endl;
    dec_bit_rt_one(K1, 1);
	std::cout << ">>: " << dec_to_string(k_sgn, K1) << std::endl;
    dec_lt_move(K1, 63);
	std::cout << "63: " << dec_to_string(k_sgn, K1) << std::endl;
    dec_lt_move(K1, 1);
	std::cout << "64: " << dec_to_string(k_sgn, K1) << std::endl;
    dec_rt_move(K1, 63);
	std::cout << "1: " << dec_to_string(k_sgn, K1) << std::endl;

    auto start = NEUNET_CHRONO_TIME_POINT;
	// for(uint64_t l = 0; l < 1000000000000000000; l++){
	//     for(uint64_t k = 0; k < 1000000000000000000; k++){
	//         for(uint64_t j = 0; j < 1000000000000000000; j++){
	//             for(uint64_t i = 0; i < 1000000000000000000; i++){
	//                 for(uint64_t h = 0; h < 1000000000000000000; h++){
	//                     for(uint64_t g = 0; g < 1000000000000000000; g++){
	//                         for(uint64_t f = 0; f < 1000000000000000000; f++){
	//                             for(uint64_t e = 0; e < 1000000000000000000; e++){
	//                                 for(uint64_t d = 0; d < 1000000000000000000; d++){
	//                                     for(uint64_t c = 0; c < 1000000000000000000; c++){
	//                                         for(uint64_t b = 0; b < 1000000000000000000; b++){
	//                                             for(uint64_t a = 0; a < 1000000000000000000; a++){
    //                                                 k03 <<= 4;
    //                                                 k03 >>= 4;
    //                                             }
    //                                         }
    //                                     }
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
	auto end = NEUNET_CHRONO_TIME_POINT;
	std::cout << "Time1: " << end - start << std::endl;
    start = NEUNET_CHRONO_TIME_POINT;
	for(uint64_t a = 0; a < 100000000; a++){
        dec_lt_move(k2, 2);
        dec_rt_move(k2, 2);
    }
	end = NEUNET_CHRONO_TIME_POINT;
	std::cout << "Time: " << end - start << std::endl;
    start = NEUNET_CHRONO_TIME_POINT;
	for(uint64_t a = 0; a < 100000000; a++){
        dec_lt_move(k1, 2);
        dec_rt_move(k1, 2);
    }
	end = NEUNET_CHRONO_TIME_POINT;
	std::cout << "Time0: " << end - start << std::endl;
	// std::cout << k1.it << std::endl;
	//std::cout << "Tell: " << dec_to_string(k_sgn, k1 / k0) << std::endl;
    dec_bit_not(K);
	std::cout << k1.it << std::endl;
	std::cout << k2.it << std::endl;
	std::cout << dec_to_string(k_sgn, K) << std::endl;
    dec_lt_move(K2, 123);
    // for(int a = 120; a < 128; a++){
    //     dec_bit_lt_one(K2);
    // 	std::cout << dec_to_string(k_sgn, K2) << std::endl;
    // 	// std::cout << K2.it << std::endl;
    // }
	std::cout << dec_to_string(k_sgn, k1) << std::endl;
	std::cout << dec_to_string(k_sgn, K2) << std::endl;
    // dec_bit_not(k3);
    auto k001 = dec_bit_and(k3, k4);
	std::cout << "3: " << dec_to_string(k_sgn, k3) << std::endl;
    // dec_bit_not(k4, true);
	std::cout << "4: " << dec_to_string(k_sgn, k4) << std::endl;
	std::cout << "k_ans: " << dec_to_string(k_sgn, k001) << std::endl;
	std::cout << "K: " << dec_to_string(k_sgn, dec_bit_val(k001)) << std::endl;
    // 10111:23 101001:41 100000010 100000001 11111110
}

