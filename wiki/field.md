# 字段

## `modulus_mode`

取模模式。默认值是 `false` ，当值是 `true` 时将从取余模式变更为取模。这个字段只对整数有效。

## `default_infinite_precision`

缺省无限精度，静态字段，默认值是 `64` 。

## `binary_bit_set`

位运算二进制比特位的个数，默认值是 `0` 。当这个字段的值是 `0` 时，位运算将没有二进制位的限制，否则计算将限定在这个值的二进制位的范围内。这个字段只对整数有效。

## `infinite_precision`

无限近似值精度，默认值是 [`default_infinite_precision`](#default_infinite_precision)。计算结果是无理数或打印数值是无限小数，那么将保留 `infinite_precision` 位小数。

## `absolute`

当前数的绝对值。属性。

## `float_point_format`

当前数的 `long double` 类型数值，当数值有效数字个数超过 `long double` 的有效数字限制时，会输出截断错误值。属性。

## `integer_format`

当前数的 `int64_t` 类型数值，当数值有效数字个数超过 `int64_t` 的有效数字限制时，会输出截断错误值。属性。
