# 函数

## 构造函数

构造函数包括直接构造，空构造以及复制移动构造。

```C++
// 空构造
net_decimal blank_num_0, blank_num_1 {};
// 直接构造
auto        num_2_321 = 2.321_d; // literal 字面量
net_decimal num_3_141 = "3.141"; // 字符串数字
net_decimal num_4_668 = 4.668;
// 复制构造
auto num_cpy_2_321 = num_2_321;
// 移动构造
auto num_mov_4_668 = std::move(num_4_668); // num_4_668 将会析构为 0
```

## `pi`

指定保留小数位，获取圆周率 $\pi$ 的近似值。

参数|IO|描述|默认值
-|-|-|-
`places`|*In*|保留 $\pi$ 的小数位数。|[`default_infinite_precision`](field.md/#default_infinite_precision)

返回 $\pi$ 的近似值

## `reduct`

约分，将分式有理数进行简化。这个函数一般用作计算提速。

## `int_part`

获取当前数的整数与浮点部分。

参数|IO|描述
-|-|-
`float_part`|*Out*|当前数的浮点部分。

返回当前数的整数部分。

获取 $2.25$ 的整数和小数部分。

```C++
using namespace std;
auto d25 = 0_d,
     i2  = (2.25_d).float_part(d25);
cout << i2 << endl; // 打印 2.25 整数部分
cout << d25 << endl; // 打印 2.25 小数部分
```

## `log`

获取当前值以 $e$ 为底的对数。

返回函数 $\ln{x}$ 的值，$x$ 为当前值。

计算 $ln{2}$ 的值并打印。

```C++
using namespace std;
cout << log(2_d) << endl;
```

## `exp`

获取当前值以 $e$ 为底的指数。

返回函数 $\exp{x}$ 的值，即 $e^x$, $x$ 为当前值。

计算 $e^{2.16}$ 的值并打印。

```C++
using namespace std;
cout << exp(2.16_d) << endl;
```

## `sin`

获取当前值的正弦值。

返回函数 $\sin{x}$ 的值，$x$ 为当前值。

计算 $\sin{1.68}$ 的值并打印

```C++
using namespace std;
cout << sin(1.68_d) << endl;
```

## `cos`

获取当前值的余弦值。

返回函数 $\cos{x}$ 的值，$x$ 为当前值。

计算 $\cos{2.72}$ 的值并打印

```C++
using namespace std;
cout << cos(2.72_d) << endl;
```

## `pow`

计算数值的幂。

参数|IO|描述
-|-|-
`_Xx`|*In*|底数
`_Yx`|*In*|次数

返回函数 $x^y$ 的值，$x$ 为底数，$y$ 为次数。

计算 $(-32) ^ {0.2}$ 的值并打印

```C++
using namespace std;
cout << pow(-32_d, 0.2_d) << endl;
```

*关于打印数值的精度更多信息，请参阅* [`infinite_precision`](field.md/#infinite_precision) *字段。*
