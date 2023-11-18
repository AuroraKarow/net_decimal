# Field

## `modulus_mode`

Modulo operation mode. Default value is `false` , while remainder mode would be changed to modulo mode for `true` value. This field is valid in integer field only.

## `default_infinite_precision`

Default precision value of infinite decimal. Static field. Default value is `64` .

## `binary_bit_set`

Binary bit count for bitwise operation. Default value is `0` . There would be no bit limit for bitwise calculation for the `0` value otherwise limited below this value. This field is valid in integer field only.

## `infinite_precision`

Precision of infinite decimal approximation. Default value is [`default_infinite_precision`](#default_infinite_precision). The answer would keep `infinite_precision` decimal places for irrational result or printing infinite decimal.

## `absolute`

Absolute value of current number. Property.

## `float_point_format`

The `long double` type of current number. It would meet truncating error when the significant figure of current number surpasses the `long double` type. Property.

## `integer_format`

The `int64_t` type of current number. It would meet truncating error when the significant figure of current number surpasses the `int64_t` type. Property.
