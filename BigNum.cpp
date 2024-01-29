//
// Created by Семён on 25.01.2024.
//

#include "BigNum.h"

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
    block_size_ = 4;
    base_ = static_cast<uint64_t>(std::pow(10, block_size_));
//    base_ = 10000u;
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
    if (exp_ > 0 && res.length() < exp_) {
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
    base_ = other.base_;
    exp_ = other.exp_;
    block_size_ = other.block_size_;
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
        accum = first.blocks_[ind] / first.base_;
        first.blocks_[ind] %= first.base_;
    }

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
                first.blocks_[ind_f] = first.base_ - 1;
            }
            first.blocks_[ind] += first.base_ - second.blocks_[ind];
        }
    }

    if (is_change_sign) {
        first.sign_ = -1;
    }

    return first;
}

BigNum::BigNum(const BigNum &other) {
    blocks_ = other.blocks_;
    sign_ = other.sign_;
    base_ = other.base_;
    exp_ = other.exp_;
    block_size_ = other.block_size_;
}

bool operator<(const BigNum &first, const BigNum &second) {
    if (first.sign_ < second.sign_) return true;

    auto str1 = first.toString(false);
    auto str2 = second.toString(false);

    for (size_t ind = 0; ind < std::max(str1.length(), str2.length()); ++ind) {
        if (ind < str1.length() && ind < str2.length()) {
            if (str1[ind] < str2[ind]) {
                return true;
            } else {
                return false;
            }
        } else if (ind < str1.length()) {
            if (str1[ind] != 0) {
                return false;
            }
        } else {
            if (str2[ind] != 0) {
                return true;
            }
        }
    }

    return false;
}

BigNum::BigNum(BigNum &&other) noexcept {
    blocks_ = std::move(other.blocks_);
    sign_ = other.sign_;
    exp_ = other.exp_;
    base_ = other.base_;
    block_size_ = other.block_size_;
}

BigNum operator*(const BigNum& first, const BigNum& second) {
    BigNum res = 0;
    for (auto block : second.blocks_) {
        uint64_t accum = 0;
        BigNum inter_term;
        for (int second_ind = 0; second_ind < first.blocks_.size(); ++second_ind) {
            auto temp_val = block * first.blocks_[second_ind] + accum;
            inter_term.blocks_.push_back(temp_val % first.base_);
//            inter_term = inter_term + (temp_val % first.base_);
            accum = (temp_val) / first.base_;
        }
        if (accum != 0) {
            inter_term.blocks_.push_back(accum);
        }
        auto str1 = inter_term.toString();
        res = res + inter_term;
        auto str2 = res.toString();
    }

    res.sign_ = first.sign_ * second.sign_;
    res.exp_ = first.exp_ + second.exp_;
    return res;
}

BigNum operator/(BigNum first, BigNum second) {
    return BigNum();
}

void swap(BigNum &lhs, BigNum &rhs) {
    std::swap(lhs.blocks_, rhs.blocks_);
    std::swap(lhs.block_size_, rhs.block_size_);
    std::swap(lhs.exp_, rhs.exp_);
    std::swap(lhs.sign_, rhs.sign_);
    std::swap(lhs.base_, rhs.base_);
}
