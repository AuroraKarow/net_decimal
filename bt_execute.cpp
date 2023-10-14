#include "iostream"
#include "net_decimal"
#include "bitset"
#include "math.h"
#include "chrono"

using namespace neunet;
#define NEUNET_CHRONO_TIME_POINT    std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()).time_since_epoch().count()


void dec_bit_ltn(net_decimal_data &D, uint64_t bit = 1, bool sta = false){
    if (D.ft.length || !D.it.length) return;
    uint64_t b = bit % NEUNET_DEC_DIG_MAX,
             B = bit / NEUNET_DEC_DIG_MAX,
             l = D.it.length,
             T = 0,
             t;
    if (!sta && (B || D.it[l - 1] >= (NEUNET_DEC_SEG_MAX >> b))) D.it.init(l + (bit / 63) + 1);
    if (D.it[l - 1] >= (NEUNET_DEC_SEG_MAX >> b)) l++;
    // std::cout << "S" << std::endl;
    if (b){
        uint64_t bas = NEUNET_DEC_SEG_MAX >> b;
        for(uint64_t m = 0; m < l; m++) {
            t = D.it[m] / bas;
            D.it[m] %= bas;
            D.it[m] <<= b;
            D.it[m] += T;
            T = t;
        }
    }
    // std::cout << "A" << std::endl;
    uint64_t Bas = NEUNET_DEC_SEG_MAX >> NEUNET_DEC_DIG_MAX;
    while (B--){
        T = 0;
        for(uint64_t m = 0; m < D.it.length; m++) {
            t = D.it[m] / Bas;
            D.it[m] %= Bas;
            D.it[m] <<= NEUNET_DEC_DIG_MAX;
            D.it[m] += T;
            T = t;
        }
    }
    l = D.it.length;
    // std::cout << "B" << std::endl;
    if (!sta) {
        while(!D.it[--l]) continue;
        D.it = D.it.sub_set(0, l);
    }
}

void dec_bit_rtn(net_decimal_data &D, uint64_t bit = 1, bool sta = false){
    if (D.ft.length || !D.it.length) return;
    uint64_t b = bit % NEUNET_DEC_DIG_MAX,
             B = bit / NEUNET_DEC_DIG_MAX,
             l = D.it.length,
             T = 0,
             t;
    // std::cout << "S" << std::endl;
    if (b){
        uint64_t bas = NEUNET_DEC_SEG_MAX >> b,
                 div = 1 << b;
        for(uint64_t m = l; m > 0; m--) {
            t = (D.it[m - 1] % div) * bas;
            D.it[m - 1] >>= b;
            D.it[m - 1] += T;
            T = t;
        }
    }
    // std::cout << "A" << std::endl;
    uint64_t Bas = NEUNET_DEC_SEG_MAX >> NEUNET_DEC_DIG_MAX,
             Div = 1 << NEUNET_DEC_DIG_MAX;
    while (B--){
        T = 0;
        for(uint64_t m = l; m > 0; m--) {
            t = (D.it[m - 1] % Div) * Bas;
            D.it[m - 1] >>= NEUNET_DEC_DIG_MAX;
            D.it[m - 1] += T;
            T = t;
        }
    }
    // std::cout << "B" << std::endl;
    if (!sta) {
        while(!D.it[--l]) continue;
        D.it = D.it.sub_set(0, l);
    }
}

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
    //if (!(D.it.length - 1)) D.it[0] >>= 1; return;
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

uint64_t rev0(uint64_t &d){
    uint64_t dn = 0;
    for(uint64_t k = 0; k < 64; k++){ 
        uint64_t a = 1, b = 1;
        if(d & (a << k)){ dn += b << (63 - k); } 
    }
    d = dn;
    return d;
}

void dec_bit_lt_ones(net_decimal_data &D, const uint64_t &b){
    if (D.ft.length || !D.it.length) return;
    int t = 0;
    for(uint64_t m = 0; m < b; m++) if(D.it[m] >= NEUNET_DEC_BIT_BAS) {
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

void dec_bit_rt_ones(net_decimal_data &D, const uint64_t &b){
    if (D.ft.length || !D.it.length) return;
    uint64_t t = 0;
    for(int64_t m = b + 1; m >= 0; m--) if(D.it[m] % 2){
        D.it[m] >>= 1;
        D.it[m] += t;
        t = NEUNET_DEC_BIT_BAS;
    } else{
        D.it[m] >>= 1;
        D.it[m] += t;
        t = 0;
    }
}

void dec_lt_move(net_decimal_data &D, int64_t B){ while(B--){ dec_bit_lt_one(D); } }

void dec_rt_move(net_decimal_data &D, int64_t B){ while(B--){ dec_bit_rt_one(D); } }

void dec_mlt(net_decimal_data &D, uint64_t B){
    uint64_t a = 0, b = dec_dig_cnt(D, true);
    D.it.init((b + B - 1) / NEUNET_DEC_DIG_MAX + 1);
    std::cout << D.it << std::endl;
    while(a++ < B){ dec_bit_lt_ones(D, (a + b - 1) / NEUNET_DEC_DIG_MAX + 1); } 
}

void dec_mrt(net_decimal_data &D, uint64_t B){
    uint64_t a = 0, b = dec_dig_cnt(D, true);
    // std::cout << D.it << std::endl;
    // std::cout << (b - B) / NEUNET_DEC_DIG_MAX << std::endl;
    while(a++ < B){ dec_bit_rt_ones(D, (b - a - 1) / NEUNET_DEC_DIG_MAX); }
    // std::cout << D.it << std::endl;
    // std::cout << (b - B) / NEUNET_DEC_DIG_MAX << std::endl;
    D.it.sub_set(0, (b - B - 1) / NEUNET_DEC_DIG_MAX);
    // std::cout << D.it << std::endl;
    // std::cout << std::endl;
}


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
    dec_bit_not(d, false, bit - 2);
    return d;
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


int main(){
    bool k_sgn = true;
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
        k1 = dec_init(k_sgn, "81234567890111213141516171819202122232"),
        k6 = dec_init(k_sgn, "1234567890111213141516171819202122232"),
        k2 = dec_init(k_sgn, "4"),
        k3 = dec_init(k_sgn, "206"),//11001110 / 306 / 100110010
        k4 = dec_init(k_sgn, "61"),//00111101 / 451 / 111000011
        k5 = dec_init(k_sgn, "0");//100000010 / 258 / 11111110

    // auto start = NEUNET_CHRONO_TIME_POINT;
    // auto ds_rev = rev(64);
	// for(int a = 0;a < 99999;a++){ dec_rev_exec0(k0, ds_rev); }
	// auto end = NEUNET_CHRONO_TIME_POINT;
	// std::cout << "Time0: " << end - start << std::endl;
	// std::cout << dec_to_string(k_sgn, k0) << std::endl;
    auto K = k1, K1 = dec_bit_k1(), K2 = dec_bit_k1();
    dec_bit_lt_one(K1);
	std::cout << "<<: " << K1 << std::endl;
    dec_bit_rt_one(K1, 1);
	std::cout << ">>: " << K1 << std::endl;
    dec_lt_move(K1, 63);
	std::cout << "63: " << K1 << std::endl;
    dec_lt_move(K1, 1);
	std::cout << "64: " << K1 << std::endl;
    dec_rt_move(K1, 63);
	std::cout << "1: " << K1 << std::endl;

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
	// std::cout << "Time1: " << end - start << std::endl;
    // dec_lt_move(k1, 2);
    start = NEUNET_CHRONO_TIME_POINT;
	for(uint64_t a = 0; a < 10000000; a++){
        dec_lt_move(k6, 2);
        dec_rt_move(k6, 2);
    }
	end = NEUNET_CHRONO_TIME_POINT;
	std::cout << "Time1: " << end - start << std::endl;
    std::cout << k6 << std::endl;

    start = NEUNET_CHRONO_TIME_POINT;
	for(uint64_t a = 0; a < 10000000; a++){
        dec_bit_ltn(k6, 2);
        // std::cout << k6 << std::endl;
        dec_bit_rtn(k6, 2);
        // std::cout << k6 << std::endl;
    }
	end = NEUNET_CHRONO_TIME_POINT;
	std::cout << "Time2: " << end - start << std::endl;
    std::cout << k6 << std::endl;

    start = NEUNET_CHRONO_TIME_POINT;
	for(uint64_t a = 0; a < 1000000; a++){
        dec_bit_ltn(k1, 100);
        // std::cout << k1 << std::endl;
        dec_bit_rtn(k1, 100);
        // dec_mrt(k1, 32);
        // std::cout << k1 << std::endl;
    }
	end = NEUNET_CHRONO_TIME_POINT;
    dec_bit_ltn(k1, 100);
	std::cout << "Time01: " << end - start << std::endl;
    std::cout << k1 << std::endl;
    // dec_mrt(k1, 32);
    dec_bit_rtn(k1, 100);

    dec_lt_move(K2, 50);
    start = NEUNET_CHRONO_TIME_POINT;
	for(uint64_t a = 0; a < 1000000; a++){
        k5 = dec_mul(k1, K2);
        k5 = dec_mul(k5, K2);
    }
	end = NEUNET_CHRONO_TIME_POINT;
	std::cout << "Time: " << end - start << std::endl;
    std::cout << k5 << std::endl;

    start = NEUNET_CHRONO_TIME_POINT;
	for(uint64_t a = 0; a < 1000000; a++){
        dec_lt_move(k1, 100);
        dec_rt_move(k1, 100);
    }
	end = NEUNET_CHRONO_TIME_POINT;
    dec_lt_move(k1, 100);
	std::cout << "Time00: " << end - start << std::endl;
    std::cout << k1 << std::endl;
    dec_rt_move(k1, 100);

	// std::cout << k1.it << std::endl;
	//std::cout << "Tell: " << dec_to_string(k_sgn, k1 / k0) << std::endl;
    dec_bit_not(K);
	std::cout << k1.it << std::endl;
	std::cout << k2.it << std::endl;
	std::cout << K << std::endl;
    dec_lt_move(K2, 123);
    // for(int a = 120; a < 128; a++){
    //     dec_bit_lt_one(K2);
    // 	std::cout << dec_to_string(k_sgn, K2) << std::endl;
    // 	// std::cout << K2.it << std::endl;
    // }
	std::cout << k3 << std::endl;
	std::cout << k4 << std::endl;
    // dec_bit_not(k3);
    auto k001 = dec_bit_and(k3, k4, true, true);
	std::cout << "3: " << k3 << std::endl;
    // dec_bit_not(k4, true);
	std::cout << "4: " << k4 << std::endl;
	std::cout << "k_ans: " << k001 << std::endl;
	std::cout << "K: " << dec_bit_val(k001) << std::endl;
    // 10111:23 101001:41 100000010 100000001 11111110
    return 0;
}

