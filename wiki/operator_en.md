# Operator

## Explicit type casting

Convert type `net_decimal` to arithmetic type.

```C++
auto i4d121  = 4.121_d;
auto i4d121d = double(i4d121);
auto i4      = int(i4d121);
```

## Assignment

Assign type of `net_decimal` or arithmetic value to `net_decimal` variable.

```C++
net_decimal dest = 0;
int i4 = 4;
dest = i4; // here dest = 4
dest = 32.116; // here dest = 32.116
dest = "163.21"; // here dest = 163.21
```

## Comparison

Compare type of `net_decimal` or arithmetic value to `net_decimal` variable.

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

## Addition

`+` `++` `+=`

Define possitive number and calculate addition and increment.

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

## Subtraction

`-` `--` `-=`

Define negative number, calculate subtraction and decrement.

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

## Multiplication

`*` `*=`

Multiplication calculation.

```C++
using namespace std;
auto i33d142 = 33.142_d,
     md92    = -.92_d;
cout << i33d142 * md92 << endl;
md92 *= i33d142;
md92 *= 2.6;
cout << md92 << endl;
```

## Division

`/` `/=`

Division calculation

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

## Remainder \& Modulo

`%` `%=`

Operator for remainder or modulo, operands are integer.

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

*For more infomation about the modulus mode, please refer to field* [`modulus_mode`](field.md/#modulus_mode) *.*

## Left shift

`<<` `<<=`

Bitwise operation for left shifting, operands are integer.

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

## Right shift

`>>` `>>=`

Bitwise operation for right shifting, operands are integer.

```C++
using namespace std;
auto itest = "192831284698325983253481263856325632"_d,
     ibit  = 12_d;
cout << (itest >> ibit) << endl;
itest >>= ibit;
itest >>= 2;
cout << (itest >> 8) << endl;
```

## Bitwise AND

`&` `&=`

Bitwise AND operation, operands are integer.

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

## Bitwise OR

`|` `|=`

Bitwise OR operation, operands are integer.

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

## Bitwise NOT

`~`

Bitwise NOT operation, operand is integer.

```C++
auto i82132 = 82132_d;
~i82132;
cout << i82132 << endl;
```

## Bitwise XOR

`^` `^=`

Bitwise XOR operation, operands are integer.

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
