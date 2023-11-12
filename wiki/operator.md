# 运算符

## 强制类型转换

将 `net_decimal` 强制转换为数值类型。

```C++
auto i4d121  = 4.121_d;
auto i4d121d = double(i4d121);
auto i4      = int(i4d121);
```

## 赋值

将 `net_decimal` 或数值类型的值，赋给 `net_decimal` 变量。

```C++
net_decimal dest = 0;
int i4 = 4;
dest = i4; // 此时 dest = 4
dest = 32.116; // 此时 dest = 32.116
dest = "163.21"; // 此时 dest = 163.21
```

## 比较运算符

将 `net_decimal` 变量与数值变量或 `net_decimal` 变量进行比较。

```C++
using namespace std;
net_decimal i4    = 4;
double      i1d41 = 1.41;
int         i6    = 6;
string      si4   = "4";
if (i4 > i1d41) cout << i4 << " is greater than " << i1d41 << endl;
if (i4 < i6) cout << i4 << " is less than " << i6 << endl;
if (i4 == si4) cout << i4 << " is equal to " << si14 << endl;
```

## 加运算符

`+` `++` `+=`

可以定义正数，计算加法以及自增。

```C++
using namespace std;
auto ai24d36 = +24.36_d,
     i10d22  = 10.22_d;
cout << ai24d36 + i10d22 << endl;
cout << ++ai24d36 << endl;
cout << i10d22++ << endl;

ai24d36 += i10d22;

double i12d4 = 12.4;
cout << ai24d36 + i12d4 << endl;
```

## 减运算符

`-` `--` `-=`

可以定义负数，计算减法以及自减。

```C++
using namespace std;
auto mi24d36 = -24.36_d,
     i10d22  = 10.22_d;
cout << i10d22 - mi24d36 << endl;
cout << --i10d22 << endl;
cout << i10d22-- << endl;

i10d22 -= mi24d36;

double i12d4 = 12.4;
cout << i10d22 + i12d4 << endl;
```

## 乘法运算符

`*` `*=`

计算乘法。

```C++
using namespace std;
auto i33d142 = 33.142_d,
     md92    = -.92_d;
cout << i33d142 * md92 << endl;
md92 *= i33d142;
md92 *= 2.6;
cout << md92 << endl;
```

## 除法运算符

`/` `/=`

计算除法。

```C++
using namespace std;
auto i76d213 = 76.213_d,
     i144    = 144_d;
auto d1453   = .1453;
cout << i76d213 / d1453 << endl;
cout << i144 / i76d213 << endl;
d1453   /= i76d213;
i76d213 /= i144;
i144    /= d1453;
cout << d1453 << endl;
cout << i76d213 << endl;
cout << i144 << endl;
```

## 余模运算符

`%` `%=`

取余或取模运算，操作对象为整数。

```C++
using namespace std;
auto mi17 = -17_d,
     i7   = 7_d;
cout << mi17 % i7 << endl;
i7.modulus_mode = true;
cout << mi17 % i7 << endl;
mi17 %= i7;
i7  %= 3;
```

*关于取模模式更多信息，请参阅* [`modulus_mode`](field.md/#modulus_mode) *字段。*

## 左移运算符

`<<` `<<=`

左移位运算，操作对象为整数。

```C++
using namespace std;
auto lsh_i22 = 22_d;
auto bit_i14 = 14_d;
cout << (lsh_i22 << bit_i14) << endl;
cout << (lsh_i22 << 4) << endl;
auto i32 = 32;
i32     <<= bit_i14;
bit_i14 <<= i32;
bit_i14 <<= lsh_i22;
```

## 右移运算符

`>>` `>>=`

右移位运算，操作对象为整数。

```C++
using namespace std;
auto itest = "192831284698325983253481263856325632"_d,
     ibit  = 12_d;
cout << (itest >> ibit) << endl;
itest >>= ibit;
itest >>= 2;
cout << (itest >> 8) << endl;
```

## 与运算符

`&` `&=`

按位与运算，操作对象为整数。

```C++
using namespace std;
auto i47 = 47_d,
     i32 = 32_d;
cout << (i47 & 16) << endl;
cout << (i47 & i32) << endl;
i32 &= 16;
i47 &= i32;
auto mask = 64;
mask &= i47;
cout << (mask & i32) << endl;
```

## 或运算符

`|` `|=`

按位或运算，操作对象为整数。

```C++
using namespace std;
auto i47 = 47_d,
     i32 = 32_d;
cout << (i47 | 16) << endl;
cout << (i47 | i32) << endl;
i32 |= 16;
i47 |= i32;
auto mask = 64;
mask |= i47;
cout << (mask | i32) << endl;
```

## 非运算符

`~`

按位取反，操作对象为整数。

```C++
auto i82132 = 82132_d;
~i82132;
cout << i82132 << endl;
```

## 异或运算符

`^` `^=`

按位异或，操作对象为整数。

```C++
using namespace std;
auto i47 = 47_d,
     i32 = 32_d;
cout << (i47 ^ 16) << endl;
cout << (i47 ^ i32) << endl;
i32 ^= 16;
i47 ^= i32;
auto mask = 64;
mask ^= i47;
cout << (mask ^ i32) << endl;
```
