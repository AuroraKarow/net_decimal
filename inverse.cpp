#include "bitset"
#include "iostream"
#include "math.h"
#include "chrono"

#define NEUNET_CHRONO_TIME_POINT    std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()).time_since_epoch().count()

int bt(int a){
	int c = 0;
	while(a > 0){
		a /= 2;
		c++;
	}
	return c;
}

uint64_t *rev(uint64_t B){
	uint64_t d = 1<<B;
	uint64_t *box = new uint64_t[d];
	uint64_t c = 1;
	*box = 0;
	while(c < d){
		//std::cout << "c: " << c << " d: " << d << std::endl;
		for(int i = 0;i < c;i++){ *(box + i + c) = *(box + i) + d / c / 2; }
		c *= 2;
		//c<<1;
	}
	return box;
}

#define bit_cnt 20

int main() {
	uint64_t d 	   = bit_cnt;
	auto bgn_tm_pt = NEUNET_CHRONO_TIME_POINT;
	for (auto i = 0ull; i < 10000; ++i) {
		auto da  = rev(d);
		delete []da;
	}
	auto dur = NEUNET_CHRONO_TIME_POINT - bgn_tm_pt;
	// auto da  = rev(d);
	// std::cout << "Bytes: " << d << std::endl;
	// for(auto a = 0;a < (1<<d);a++){
	// 	std::cout << a << " - " << std::bitset<bit_cnt>(a) << " -> " << std::bitset<bit_cnt>(*(da + a)) << " - " << *(da + a) << std::endl;
	// }
	// delete []da;
	std::cout << dur << "ms" << std::endl;
	
	return 0;
}