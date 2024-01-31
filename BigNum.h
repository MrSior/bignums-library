#ifndef BIGNUMSLIBRARY_BIGNUM_H
#define BIGNUMSLIBRARY_BIGNUM_H

#include <vector>
#include <bitset>
#include <iostream>
#include <type_traits>
#include <concepts>
#include "BigNumError.h"

template <typename T>
concept Number = std::integral<T> || std::floating_point<T>;

template <typename T>
concept ConvertableToString = std::constructible_from<std::string, T>;

template<typename T>
concept HasToString = requires (T t) {
    { t.toString() } -> std::same_as<std::string>;
};

class BigNum {
private:
    std::vector<uint16_t> blocks_;
    int sign_;
    static const uint64_t base_ = 10000;
    int32_t exp_;
    int block_size_;

    void Init();
    void Init(std::string str);


    void SeparateToBlocks(const std::string& str);
    void RemoveInsignificantZeroes();
public:
    BigNum();

    template<HasToString T>
    BigNum(T val) {
        Init(val.toString());
    }

    template<ConvertableToString T>
    BigNum(T val) {
        Init(std::string(val));
    }

    template<Number T>
    BigNum(T val) {
        Init(std::to_string(val));
    }

    BigNum(const BigNum& other);
    BigNum(BigNum&& other) noexcept;

    BigNum& operator=(BigNum&& other) noexcept;
    BigNum& operator=(const BigNum& other);

    void printBlocks();

    void ShiftToExp(int32_t new_exp);

    [[nodiscard]] std::string toString(bool with_sign = true, bool with_dote = true) const;
    int getSign() { return sign_; }
    std::vector<uint16_t> getBlocks() { return blocks_; }

    friend void swap(BigNum& lhs, BigNum& rhs);

    friend BigNum operator+(BigNum, BigNum);
    friend BigNum operator-(BigNum, BigNum);
    friend BigNum operator/(BigNum, const BigNum&);
    friend BigNum operator*(const BigNum&, const BigNum&);

    friend bool operator<(const BigNum&, const BigNum&);
    friend bool operator<=(const BigNum&, const BigNum&);
    friend bool operator!=(const BigNum&, const BigNum&);
    friend bool operator==(const BigNum&, const BigNum&);
    friend bool operator>(const BigNum&, const BigNum&);
    friend bool operator>=(const BigNum&, const BigNum&);

};

static BigNum operator"" _bn(const char* val) {
    return {val};
}

#endif //BIGNUMSLIBRARY_BIGNUM_H
