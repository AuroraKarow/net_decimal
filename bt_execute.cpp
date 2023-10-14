#include "iostream"
#include "net_decimal"
#include "bitset"
#include "math.h"
#include "chrono"

using namespace neunet;
#define NEUNET_CHRONO_TIME_POINT    std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()).time_since_epoch().count()


void dec_lt_move_one(net_decimal_data &D){
    bool k_sgn = false;
    uint64_t k1 = 5000000000000000000, k = 1;
    k <<= 63;
    if(!dec_dig_cnt(D, false)){
        if(D.it[D.it.length - 1] >= k){
            D.it[D.it.length - 1] -= k;
            D.it.init(D.it.length + 1);
            D.it[D.it.length - 1] += k;
            int t = 0;
            for(uint64_t m = 0; m < D.it.length - 1; m++){
                if(D.it[m] >= k1){ 
                    D.it[m] -= k1;
                    D.it[m] <<= 1;
                    D.it[m] += t;
                    t = 1;
                }
                else{
                    D.it[m] <<= 1;
                    D.it[m] += t;
                    t = 0;
                }
            }
        }
        else{
            int t = 0;
            for(uint64_t m = 0; m < D.it.length; m++){
                if(D.it[m] >= k1){ 
                    D.it[m] -= k1;
                    D.it[m] <<= 1;
                    D.it[m] += t;
                    t = 1;
                }
                else{
                    D.it[m] <<= 1;
                    D.it[m] += t;
                    t = 0;
                }
            }
            if(t != 0){
                D.it.init(D.it.length + 1);
                D.it[D.it.length - 1]++;  
            }
        }
        // std::cout << "dec_lt" << std::endl;
    }
}

void dec_rt_move_one(net_decimal_data &D){
    bool k_sgn = false;
    uint64_t k0 = 0, k1 = 5000000000000000000, k = 1;
    k <<= 63;
    if(!dec_dig_cnt(D, false)){
        if(D.it[D.it.length - 1] >= k){
            D.it[D.it.length - 1] -= k;
            uint64_t t = 0;
            for(uint64_t m = D.it.length; m > 0; m--){
                if(D.it[m - 1] % 2 == 1){ 
                    D.it[m - 1] >>= 1;
                    D.it[m - 1] += t;
                    t = k1;
                }
                else{
                    D.it[m - 1] >>= 1;
                    D.it[m - 1] += t;
                    t = 0;
                }
            }
            if(D.it[D.it.length - 1] > 0){ D.it.init(D.it.length + 1); }
            D.it[D.it.length - 1] += k;
        }
        else{
            uint64_t t = 0;
            for(uint64_t m = D.it.length; m > 0; m--){
                if(D.it[m - 1] % 2 == 1){ 
                    D.it[m - 1] >>= 1;
                    D.it[m - 1] += t;
                    t = k1;
                }
                else{
                    D.it[m - 1] >>= 1;
                    D.it[m - 1] += t;
                    t = 0;
                }
            }
            if(!D.it[D.it.length - 1]){ D.it.init(D.it.length - 1, true); }
        }
        // std::cout << "dec_rt" << std::endl;
    }
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

void dec_lt_move(net_decimal_data &D, int64_t B){ while(B--){ dec_lt_move_one(D); } }

void dec_rt_move(net_decimal_data &D, int64_t B){ while(B--){ dec_rt_move_one(D); } }

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
            dec_lt_move_one(kp);
            dec_rt_move_one(D1);
            dec_rt_move_one(d1);
        }
    }
    if(type == 1){
        while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
            if(!(!(D1.it[0] % 2) && !(d1.it[0] % 2))){ k = dec_add(sgn, k, sgn, kp, sgn); }
            dec_lt_move_one(kp);
            dec_rt_move_one(D1);
            dec_rt_move_one(d1);
        }
    }
    if(type == 2){
        while(dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
            if(!(D1.it[0] % 2)){ k = dec_add(sgn, k, sgn, kp, sgn); }
            dec_lt_move_one(kp);
            dec_rt_move_one(D1);
        }
    }
    if(type == 3){
        while(dec_comp(d1, k0) != NEUNET_DEC_CMP_EQL && dec_comp(D1, k0) != NEUNET_DEC_CMP_EQL){
            if((D1.it[0] % 2 && !(d1.it[0] % 2)) || (!(D1.it[0] % 2) && d1.it[0] % 2)){ k = dec_add(sgn, k, sgn, kp, sgn); }
            dec_lt_move_one(kp);
            dec_rt_move_one(D1);
            dec_rt_move_one(d1);
        }
    }
    return k;
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


void main(){
    bool k_sgn = false;
    uint64_t k01 = 1;
    k01 <<= 63;
    uint64_t k02 = 1;
    k02 <<= 62;
    uint64_t k03 = 1;
    k03 <<= 59;
    uint64_t k = k01 + k02 + k03;
	std::cout << "k:" << std::endl;
	std::cout << k01 << std::endl;
	std::cout << k02 << std::endl;
	std::cout << k03 << std::endl;
	std::cout << k << std::endl;
    // auto k0 = dec_init(k_sgn, "190000000000000000019000000000000000001900000000000000000120000000000000000019");
    auto k0 = dec_init(k_sgn, "1"),
        k1 = dec_init(k_sgn, "10000000000000000000000000000000000000"),
        k2 = dec_init(k_sgn, "160");

    // auto start = NEUNET_CHRONO_TIME_POINT;
    // auto ds_rev = rev(64);
	// for(int a = 0;a < 99999;a++){ dec_rev_exec0(k0, ds_rev); }
	// auto end = NEUNET_CHRONO_TIME_POINT;
	// std::cout << "Time0: " << end - start << std::endl;
	// std::cout << dec_to_string(k_sgn, k0) << std::endl;

    auto start = NEUNET_CHRONO_TIME_POINT;
	for(int a = 0;a < 1000000;a++){
        dec_lt_move(k1, 2);
        dec_rt_move(k1, 2); 
    }
	auto end = NEUNET_CHRONO_TIME_POINT;
	std::cout << "Time: " << end - start << std::endl;
	// std::cout << k1.it << std::endl;
	//std::cout << "Tell: " << dec_to_string(k_sgn, k1 / k0) << std::endl;
    auto K = dec_exec_bit(2, k2, k2);
	std::cout << k1.it << std::endl;
	std::cout << k2.it << std::endl;
	std::cout << dec_to_string(k_sgn, K) << std::endl;
}

