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
    static const int block_size_ = 4;
    static const int32_t division_accuracy = 100;

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
    int32_t getExp() { return exp_; }
    void setSign(int sign) { this->sign_ = sign; }
    void setExp(int32_t exp) { this->exp_ = exp; }
    bool isOdd() const { return blocks_.front() & 1; }

    void evalf(int32_t precision);

    friend void swap(BigNum& lhs, BigNum& rhs);

    friend BigNum operator+(BigNum, BigNum);
    friend BigNum operator-(BigNum, BigNum);
    friend BigNum operator/(BigNum, const BigNum&);
    friend BigNum operator*(const BigNum&, const BigNum&);
    friend BigNum BigNumDiv(BigNum first, const BigNum& second, int32_t precision, bool debug);

//    friend BigNum Division(BigNum, const BigNum&, uint16_t precision);

    friend bool operator<(const BigNum&, const BigNum&);
    friend bool operator<=(const BigNum&, const BigNum&);
    friend bool operator!=(const BigNum&, const BigNum&);
    friend bool operator==(const BigNum&, const BigNum&);
    friend bool operator>(const BigNum&, const BigNum&);
    friend bool operator>=(const BigNum&, const BigNum&);
};

BigNum operator"" _bn(const char* val);
BigNum CalcPi(int32_t precision);

#endif //BIGNUMSLIBRARY_BIGNUM_H
