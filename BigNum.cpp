//
// Created by Семён on 25.01.2024.
//

#include "BigNum.h"

#include <utility>

BigNum::BigNum() {
    Init();
}

void BigNum::printBlocks() {
    std::cout << "Sign = " << sign_ << std::endl;
    for (auto &item: blocks_) {
        std::cout << item << " ";
    }
    std::cout << std::endl;
    std::cout << "Exp = " << exp_ << std::endl;
}

void BigNum::Init() {
    blocks_ = std::vector<uint16_t>(0);
    sign_ = 1;
    exp_ = 0;
}

void BigNum::Init(std::string str) {
    try {
        if (str.empty()) {
            throw BigNumError("empty string was given to constructor");
        }
        Init();

        if (str[0] == '-') {
            sign_ = -1;
            str.erase(str.begin());
        }
        while (str.front() == '0') {
            str.erase(str.begin());
        }
        auto dote_itr = std::find(str.begin(), str.end(), '.');
        if (dote_itr != str.end()) {
            while (str.back() == '0') {
                str.pop_back();
            }
            exp_ = std::distance(str.end(), std::find(str.begin(), str.end(), '.')) + 1;
            str.erase(dote_itr);
        } else {
            while (str.back() == '0') {
                str.pop_back();
                ++exp_;
            }
        }

        if (str.empty()) {
            str = "0";
            exp_ = 0;
        }

        std::reverse(str.begin(), str.end());

        SeparateToBlocks(str);
    } catch (const BigNumError &error) {
        std::cout << "ERROR: " << error.what() << std::endl;
        exit(1);
    }
}

std::string BigNum::toString(bool with_sign, bool with_dote) const {
    std::string res;
    for (auto itr = blocks_.end() - 1; itr >= blocks_.begin(); --itr) {
        std::string num = std::to_string(*itr);
        if (itr == blocks_.end() - 1) {
            res += num;
        } else {
            res += std::string(block_size_ - num.length(), '0') + num;
        }
    }

    if (exp_ < 0 && res.length() < -exp_) {
        res = std::string((-exp_) - res.length(), '0') + res;
    }
    if (exp_ > 0) {
        res += std::string(exp_, '0');
    }

    if (with_dote && exp_ < 0) {
        res.insert(res.begin() + (res.length() + exp_), '.');
        if (res.front() == '.') {
            res = "0" + res;
        }
    }

    if (with_sign) {
        res = (sign_ == -1 ? "-" : "") + res;
    }
    return res;
}

void BigNum::SeparateToBlocks(const std::string &str) {
    for (size_t ind = 0; ind < str.length(); ind += block_size_) {
        std::string block;
        if (str.length() - ind > block_size_) {
            block = str.substr(ind, block_size_);
        } else {
            block = str.substr(ind, str.length() - ind);
        }
        std::reverse(block.begin(), block.end());
        blocks_.push_back(std::stoi(block));
    }
}

void BigNum::ShiftToExp(int32_t new_exp) {
    if (new_exp == exp_) return;
    if (new_exp > exp_) {
        throw BigNumError("ShiftToExp can only shift to lower exponent");
    }
    decltype(exp_) diff_exp = exp_ - new_exp;

    std::string str;
    for (auto itr = blocks_.end() - 1; itr >= blocks_.begin(); --itr) {
        std::string num = std::to_string(*itr);
        if (itr == blocks_.end() - 1) {
            str += num;
        } else {
            str += std::string(block_size_ - num.length(), '0') + num;
        }

        if (itr == blocks_.begin()) {
            str += std::string(diff_exp, '0');
        }
    }

    exp_ = new_exp;
    blocks_.clear();
    std::reverse(str.begin(), str.end());
    SeparateToBlocks(str);
}

BigNum &BigNum::operator=(BigNum &&other) noexcept {
    if (this == &other) {
        return *this;
    }
    blocks_ = std::move(other.blocks_);
    sign_ = other.sign_;
    exp_ = other.exp_;
    return *this;
}

BigNum &BigNum::operator=(const BigNum &other) {
    if (this == &other) {
        return *this;
    }
    *this = BigNum(other);
    return *this;
}

BigNum operator+(BigNum first, BigNum second) {
    if (first.sign_ != second.sign_) {
        if (first.sign_ == -1) {
            first.sign_ = 1;
            return second - first;
        }
        second.sign_ = 1;
        return first - second;
    }

    first.ShiftToExp(std::min(first.exp_, second.exp_));
    second.ShiftToExp(std::min(first.exp_, second.exp_));

    size_t new_num_len = std::max(first.blocks_.size(), second.blocks_.size());
    int64_t accum = 0;
    for (size_t ind = 0; ind < new_num_len || accum > 0; ++ind) {
        if (ind == first.blocks_.size()) {
            first.blocks_.push_back(0);
        }
        if (ind < second.blocks_.size()) {
            first.blocks_[ind] += second.blocks_[ind] + accum;
        } else {
            first.blocks_[ind] += accum;
        }
        accum = first.blocks_[ind] / BigNum::base_;
        first.blocks_[ind] %= BigNum::base_;
    }

    first.RemoveInsignificantZeroes();
    return first;
}

BigNum operator-(BigNum first, BigNum second) {
    if (first.sign_ == -1 && second.sign_ == -1) {
        second.sign_ = 1;
        first.sign_ = 1;
        return second - first;
    }
    if (first.sign_ == 1 && second.sign_ == -1) {
        second.sign_ = 1;
        return first + second;
    }
    if (first.sign_ == -1 && second.sign_ == 1) {
        first.sign_ = 1;
        BigNum res = first + second;
        res.sign_ = -1;
        return res;
    }

    bool is_change_sign = false;
    if (first < second) {
         is_change_sign = true;
         BigNum tmp = first;
         first = second;
         second = tmp;
    }

    first.ShiftToExp(std::min(first.exp_, second.exp_));
    second.ShiftToExp(std::min(first.exp_, second.exp_));

    for (size_t ind = 0; ind < second.blocks_.size(); ++ind) {
        if (first.blocks_[ind] >= second.blocks_[ind]) {
            first.blocks_[ind] -= second.blocks_[ind];
        } else {
            for (size_t ind_f = ind + 1; ind_f < first.blocks_.size(); ++ind_f) {
                if (first.blocks_[ind_f] != 0) {
                    --first.blocks_[ind_f];
                    break;
                }
                first.blocks_[ind_f] = BigNum::base_ - 1;
            }
            first.blocks_[ind] += BigNum::base_ - second.blocks_[ind];
        }
    }

    if (is_change_sign) {
        first.sign_ = -1;
    }

    first.RemoveInsignificantZeroes();
    return first;
}

BigNum::BigNum(const BigNum &other) {
    blocks_ = other.blocks_;
    sign_ = other.sign_;
    exp_ = other.exp_;
}

bool operator<(const BigNum &first, const BigNum &second) {
    if (first.sign_ < second.sign_) return true;
    if (first.sign_ > second.sign_) return false;

    auto str1 = first.toString(false);
    auto str2 = second.toString(false);

    auto itr1 = std::find(str1.begin(), str1.end(), '.');
    auto itr2 = std::find(str2.begin(), str2.end(), '.');
    if (std::distance(str1.begin(), itr1) < std::distance(str2.begin(), itr2)) return true;
    if (std::distance(str1.begin(), itr1) > std::distance(str2.begin(), itr2)) return false;

    for (size_t ind = 0; str1.begin() + ind != itr1; ++ind) {
        if (str1[ind] < str2[ind]) return true;
        if (str1[ind] > str2[ind]) return false;
    }

    while (itr1 < str1.end() && itr2 < str2.end()) {
        if (*itr1 > *itr2) return false;
        if (*itr1 < *itr2) return true;
        ++itr1;
        ++itr2;
    }

    return false;
}

BigNum::BigNum(BigNum &&other) noexcept {
    blocks_ = std::move(other.blocks_);
    sign_ = other.sign_;
    exp_ = other.exp_;
}

BigNum operator*(const BigNum& first, const BigNum& second) {
    if (first < second) {
        return second * first;
    }
    BigNum res = 0;
    int cur_block_ind = 0;
    for (auto block : second.blocks_) {
        uint64_t accum = 0;
        BigNum inter_term;
        for (auto block_in_first : first.blocks_) {
            auto temp_val = block * block_in_first + accum;
            inter_term.blocks_.push_back(temp_val % BigNum::base_);
            accum = (temp_val) / BigNum::base_;
        }
        if (accum != 0) {
            inter_term.blocks_.push_back(accum);
        }
        inter_term.exp_ = cur_block_ind++ * BigNum::block_size_;
        res = res + inter_term;
    }

    res.sign_ = first.sign_ * second.sign_;
    res.exp_ = first.exp_ + second.exp_;

    res.RemoveInsignificantZeroes();
    return res;
}

BigNum Division(BigNum first, const BigNum& second, uint16_t precision) {
    return BigNumDiv(first, second, precision);
}

void swap(BigNum &lhs, BigNum &rhs) {
    std::swap(lhs.blocks_, rhs.blocks_);
    std::swap(lhs.exp_, rhs.exp_);
    std::swap(lhs.sign_, rhs.sign_);
}

bool operator!=(const BigNum& first, const BigNum& second) {
    return !(first == second);
}

bool operator==(const BigNum& first, const BigNum& second) {
    auto first_str = first.toString();
    auto second_str = second.toString();

    auto itr1 = first_str.begin();
    auto itr2 = second_str.begin();

    while (itr1 < first_str.end() && itr2 < second_str.end()) {
        if (*itr1 != *itr2) return false;
        itr1++;
        itr2++;
    }

    auto is_only_zeroes = [](decltype(itr1)& itr, decltype(first_str)& str) {
        if (itr == str.end()) return true;

        if (itr < str.end() && *itr != '.') return false;
        ++itr;
        while (itr < str.end()) {
            if (*itr != '0') return false;
            ++itr;
        }
        return true;
    };

    return is_only_zeroes(itr1, first_str) && is_only_zeroes(itr2, second_str);
}

bool operator>(const BigNum& first, const BigNum& second) {
    return second < first && second != first;
}

bool operator>=(const BigNum& first, const BigNum& second) {
    return second < first || second == first;
}

bool operator<=(const BigNum& first, const BigNum& second) {
    return first < second || first == second;
}

void BigNum::RemoveInsignificantZeroes() {
    auto str = toString();
    std::cout << "LOG(RemoveInsignificantZeroes): was = " << str;

    bool isAllZero = true;
    for (auto& elem : blocks_) {
        if (elem != 0) {
            isAllZero = false;
        }
    }
    if (isAllZero) {
        Init("0");
        str = toString();
        std::cout << "  become = " << str << std::endl;
        return;
    }

    if (exp_ >= 0) {
        if (blocks_.front() == 0) {
            blocks_.erase(blocks_.begin());
        }
        if (blocks_.empty()) {
            blocks_.push_back(0);
        }
    } else {
        size_t covered_blocks = (-exp_) / block_size_ + ((-exp_) % block_size_ == 0 ? 0 : 1);
        while (blocks_.size() > covered_blocks) {
            if (blocks_.back() != 0) break;
            blocks_.pop_back();
        }
    }

    str = toString();
    std::cout << "  become = " << str << std::endl;
}

BigNum BigNumDiv(BigNum first, const BigNum& second, int32_t cur_precision) {
    if (cur_precision < 0) return 0;
    if (first == 0) return 0;

    BigNum res = 0;
    BigNum divider = 0;
    std::cout << std::endl;
    std::cout << "LOG: first = " << first.toString() << "   second = " << second.toString() << std::endl;
    std::cout << "LOG: calculating divider ..." << std::endl;
    while ((divider + 1) * second <= first) {
        divider = divider + 1;
    }
    std::cout << "LOG: divider = " << divider.toString() << std::endl;

    if (divider == 0) {
        first = first * 10;
        std::cout << "LOG: divider == 0   ==>   res = " << first.toString() << " / " << second.toString() << std::endl;
        res = BigNumDiv(std::move(first), second, cur_precision - 1);
        --res.exp_;
    } else {
        first = first - second * divider;
        std::cout << "LOG: divider == "<< divider.toString() <<"  ==>   res = " << first.toString() << " / " << second.toString() << std::endl;
        res = divider + BigNumDiv(std::move(first), second, cur_precision);
    }

    res.RemoveInsignificantZeroes();
    return res;
}

BigNum operator/(BigNum first, const BigNum& second) {
    return BigNumDiv(std::move(first), second, BigNum::division_accuracy);
}
