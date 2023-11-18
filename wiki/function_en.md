# Function

## Constructor

Constructors include direct, blank and copy & move constructors.

```C++
// blank constructor
net_decimal blank_num_0, blank_num_1 {};
// direct constructor
auto        num_2_321 = 2.321_d; // literal value
net_decimal num_3_141 = "3.141"; // string number
net_decimal num_4_668 = 4.668;
// copy constructing
auto num_cpy_2_321 = num_2_321;
// move constructing
auto num_mov_4_668 = std::move(num_4_668); // num_4_668 will destruct to 0
```

## `pi`

Specify the decimal place to keep, get Pi $\pi$ approximation.

Parameters|IO|Description|Default value
-|-|-|-
`places`|*In*|Decimal spaces for $\pi$ to keep.|[`default_infinite_precision`](field.md/#default_infinite_precision)

Return the approximation of $\pi$ .

## `reduct`

Reduction. Simplify the rational number. This function would speed up the calculation procedure.

## `int_part`

Get integer & decimal part of current number.

Parameters|IO|Description
-|-|-
`float_part`|*Out*|Decimal part of current number

Return integer part of current number.

Get integer & decimal part of $2.25$ .

```C++
using namespace std;
auto d25 = 0_d, // decimal part
     i2  = (2.25_d).float_part(d25); // integer part
cout << i2 << endl;
cout << d25 << endl;
```

## `log`

Get logarithm of current number with base number $e$ .

Return the value of function $\ln{x}$ , $x$ is current number.

Need $ln{2}$ and print it.

```C++
using namespace std;
cout << log(2_d) << endl;
```

## `exp`

Get power of current number with base number $e$ .

Return the value of function $\exp{x}$ , the $e^x$ , $x$ is current number.

Need $e^{2.16}$ and print it.

```C++
using namespace std;
cout << exp(2.16_d) << endl;
```

## `sin`

Get sine value of current number.

Return the value of function $\sin{x}$ , $x$ is current number.

Need $\sin{1.68}$ and print it.

```C++
using namespace std;
cout << sin(1.68_d) << endl;
```

## `cos`

Get cosine value of current number.

Return the value of function $\cos{x}$ , $x$ is current number.

Need $\cos{2.72}$ and print it.

```C++
using namespace std;
cout << cos(2.72_d) << endl;
```

## `pow`

Power calculating.

Parameters|IO|Description
-|-|-
`_Xx`|*In*|Base
`_Yx`|*In*|Times

Return the value of $x^y$ , $x$ is base, $y$ is times.

Need $(-32) ^ {0.2}$ and print it.

```C++
using namespace std;
cout << pow(-32_d, 0.2_d) << endl;
```

*For more information about the precision of printing value, please refer to field* [`infinite_precision`](field.md/#infinite_precision) *.*
