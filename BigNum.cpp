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
//    if (new_exp == exp_) return;
//    if (new_exp > exp_) {
//        throw BigNumError("ShiftToExp can only shift to lower exponent");
//    }
//    decltype(exp_) diff_exp = exp_ - new_exp;
//
//    std::string str;
//    for (auto itr = blocks_.end() - 1; itr >= blocks_.begin(); --itr) {
//        std::string num = std::to_string(*itr);
//        if (itr == blocks_.end() - 1) {
//            str += num;
//        } else {
//            str += std::string(block_size_ - num.length(), '0') + num;
//        }
//
//        if (itr == blocks_.begin()) {
//            str += std::string(diff_exp, '0');
//        }
//    }
//
//    exp_ = new_exp;
//    blocks_.clear();
//    std::reverse(str.begin(), str.end());
//    SeparateToBlocks(str);

    if (new_exp == exp_) return;
    if (new_exp > exp_) {
        throw BigNumError("ShiftToExp can only shift to lower exponent");
    }
    decltype(exp_) diff_exp = exp_ - new_exp;

    while (diff_exp > 0) {
        if (diff_exp >= 4) {
            blocks_.insert(blocks_.begin(), 0);
            diff_exp -= 4;
        } else {
            int32_t mult = std::pow(10, diff_exp);
            int64_t acum = 0;
            for (auto& block : blocks_) {
                int64_t temp_val = block * mult + acum;
                acum = (temp_val) / BigNum::base_;
                block = (temp_val) % BigNum::base_;
            }
            if (acum != 0) {
                blocks_.push_back(acum);
            }
            diff_exp = 0;
        }
    }

    exp_ = new_exp;
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
//    if (first.sign_ < second.sign_) return true;
//    if (first.sign_ > second.sign_) return false;
//
//    auto str1 = first.toString(false);
//    auto str2 = second.toString(false);
//
//    auto itr1 = std::find(str1.begin(), str1.end(), '.');
//    auto itr2 = std::find(str2.begin(), str2.end(), '.');
//    if (std::distance(str1.begin(), itr1) < std::distance(str2.begin(), itr2)) return true;
//    if (std::distance(str1.begin(), itr1) > std::distance(str2.begin(), itr2)) return false;
//
//    for (size_t ind = 0; str1.begin() + ind != itr1; ++ind) {
//        if (str1[ind] < str2[ind]) return true;
//        if (str1[ind] > str2[ind]) return false;
//    }
//
//    while (itr1 < str1.end() && itr2 < str2.end()) {
//        if (*itr1 > *itr2) return false;
//        if (*itr1 < *itr2) return true;
//        ++itr1;
//        ++itr2;
//    }
//
//    return false;
    if (first.sign_ < second.sign_) return true;
    if (first.sign_ > second.sign_) return false;

    auto max_precision = std::min(first.exp_, second.exp_);
    auto lh = first;
    auto rh = second;
    lh.ShiftToExp(max_precision);
    rh.ShiftToExp(max_precision);

    auto str_lh = lh.toString();
    auto str_rh = rh.toString();
    if (str_lh.length() < str_rh.length()) return true;
    if (str_lh.length() > str_rh.length()) return false;

    for (size_t ind = 0; ind < str_lh.length(); ++ind) {
        if (str_lh[ind] < str_rh[ind]) {
            return true;
        }
        if (str_lh[ind] > str_rh[ind]) {
            return false;
        }
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
    res.exp_ += first.exp_ + second.exp_;

    res.RemoveInsignificantZeroes();
    return res;
}

//BigNum Division(BigNum first, const BigNum& second, uint16_t precision) {
//    return BigNumDiv(std::move(first), second, precision);
//}

void swap(BigNum &lhs, BigNum &rhs) {
    std::swap(lhs.blocks_, rhs.blocks_);
    std::swap(lhs.exp_, rhs.exp_);
    std::swap(lhs.sign_, rhs.sign_);
}

bool operator!=(const BigNum& first, const BigNum& second) {
    return !(first == second);
}

bool operator==(const BigNum& first, const BigNum& second) {
//    auto first_str = first.toString();
//    auto second_str = second.toString();
//
//    auto itr1 = first_str.begin();
//    auto itr2 = second_str.begin();
//
//    while (itr1 < first_str.end() && itr2 < second_str.end()) {
//        if (*itr1 != *itr2) return false;
//        itr1++;
//        itr2++;
//    }
//
//    auto is_only_zeroes = [](decltype(itr1)& itr, decltype(first_str)& str) {
//        if (itr == str.end()) return true;
//
//        if (itr < str.end() && *itr != '.') return false;
//        ++itr;
//        while (itr < str.end()) {
//            if (*itr != '0') return false;
//            ++itr;
//        }
//        return true;
//    };
//
//    return is_only_zeroes(itr1, first_str) && is_only_zeroes(itr2, second_str);
    auto max_precision = std::min(first.exp_, second.exp_);
    auto lh = first;
    auto rh = second;
    lh.ShiftToExp(max_precision);
    rh.ShiftToExp(max_precision);

    auto str_lh = lh.toString();
    auto str_rh = rh.toString();
    if (str_lh.length() != str_rh.length()) return false;

    for (size_t ind = 0; ind < str_lh.length(); ++ind) {
        if (str_lh[ind] != str_rh[ind]) {
            return false;
        }
    }
    return true;
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
    bool isAllZero = true;
    for (auto& elem : blocks_) {
        if (elem != 0) {
            isAllZero = false;
        }
    }
    if (isAllZero) {
        Init("0");
        return;
    }

    if (exp_ >= 0) {
        while (blocks_.back() == 0) {
            blocks_.pop_back();
        }
        while (blocks_.front() == 0) {
            blocks_.erase(blocks_.begin());
            exp_ += BigNum::block_size_;
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
}

BigNum BigNumDiv(BigNum first, const BigNum& second, int32_t precision, bool debug = false) {
    if (precision < 0) return 0;
    if (first == 0) return 0;

    BigNum res = 0;
    BigNum divider = 0;
    BigNum lh = 0;
    BigNum rh = first + 1;

    if (debug) {
        std::cout << std::endl;
        std::cout << "LOG: first = " << first.toString() << "   second = " << second.toString() << std::endl;
        std::cout << "LOG: calculating divider ..." << std::endl;
    }

    if (first >= second) {
        int cnt = 0;
        while (rh - lh > 1) {
            auto mid = lh + (rh - lh) * 0.5_bn;
            mid.evalf(0);

            if (debug) {
                cnt++;
                if (cnt > 300) {
                    exit(1);
                }
                std::cout << "LOG (BINSEARCH): lh = " << lh.toString() << "\tmid = " << mid.toString() << "\trh = "
                          << rh.toString() << std::endl;
            }

            if (second * mid > first) {
                rh = mid;
            } else if (second * mid <= first) {
                lh = mid;
            }
        }
        divider = lh;
    }

//    while ((divider + 1) * second <= first) {
//        std::cout << divider.toString() << std::endl;
//        if (divider.toString() == "31415926533") {
//            std::cout << ((divider + 1) * second).toString() << '\n';
//        }
//        divider = divider + 1;
//    }

    if (debug) {
        std::cout << "LOG: divider = " << divider.toString() << std::endl;
    }
    if (divider == 0) {
        first = first * 10;
        if (debug) {
            std::cout << "LOG: divider == 0   ==>   res = " << first.toString() << " / " << second.toString()
                      << std::endl;
        }
        res = BigNumDiv(std::move(first), second, precision - 1, debug);
        --res.exp_;
    } else {
        first = first - second * divider;
        if (debug) {
            std::cout << "LOG: divider == " << divider.toString() << "  ==>   res = " << first.toString() << " / "
                      << second.toString() << std::endl;
        }
        res = divider + BigNumDiv(std::move(first), second, precision, debug);
    }

    res.RemoveInsignificantZeroes();
    return res;
}

BigNum operator/(BigNum first, const BigNum& second) {
    return BigNumDiv(std::move(first), second, BigNum::division_accuracy);
}

void BigNum::evalf(int32_t precision) {
    while (precision > exp_) {
        if (precision - exp_ >= 4) {
            blocks_.erase(blocks_.begin());
            exp_ += 4;
        } else {
            *blocks_.begin() -= blocks_.front() % static_cast<int>(std::pow(10, precision - exp_));
            break;
        }
    }
}

BigNum operator"" _bn(const char* val) {
    return {val};
}

BigNum CalcPi(int32_t precision) {

//    auto arctan = [](const BigNum& x, int steps) {
//        BigNum res = x;
//        auto numer = x * x * x;
//        int denom = 3;
////        std::cout << x.toString() << std::endl;
//        for (int i = 1; i < steps; ++i, denom += 2, numer = numer * x * x) {
//            auto sum = BigNumDiv(numer, denom, 70, false);
////            std::cout << numer.toString() << " / " << denom << " = " << sum.toString() << std::endl;
//            if (i % 2 == 0) {
//                res = res + sum;
//            } else {
//                res = res - sum;
//            }
//            std::cout << "Step = " << i << "\n";
//        }
//        return res;
//    };
//
//    auto arc1 = arctan(1_bn / 5_bn, 200);
//    auto arc2 = arctan(BigNumDiv(1_bn, 239_bn, 70), 40);
//    return 4 * (4 * arc1 - arc2);
    ++precision;

    auto math_sqrt = [](const BigNum& num) {
        BigNum lh = 0;
        BigNum rh = num + 1;
        BigNum eps = 1;
        eps.setExp(-70);
        while (rh - lh > eps) {
            auto mid = lh + (rh - lh) * 0.5;

            if (mid * mid < num) {
                lh = mid;
            } else {
                rh = mid;
            }
        }
        return lh + (rh - lh) * 0.5;
    };

    auto sqrt = [&math_sqrt](const BigNum& n, const BigNum& one) {
        BigNum floating_point_precision = "10000000000000000";
        auto n_float = BigNumDiv(n * floating_point_precision, one, 0) / floating_point_precision;
        auto tmp = (floating_point_precision * math_sqrt(n_float) * one);
        tmp.evalf(0);
        auto x = BigNumDiv(tmp, floating_point_precision, 0);
        auto n_one = n * one;
        while (1) {
            auto x_old = x;
            x = BigNumDiv(x + BigNumDiv(n_one , x, 0), 2, 0);
            if (x == x_old) {
                break;
            }
        }
        return x;
    };

    BigNum C = 640320;
    BigNum C3_OVER_24 = "10939058860032000";

    struct triple {
        BigNum first;
        BigNum second;
        BigNum third;
    };

    std::function<triple(int64_t, int64_t)> bs;
    bs = [&bs, &C3_OVER_24](const int64_t& a, const int64_t& b) -> triple {
        BigNum Pab = 0;
        BigNum Qab = 0;
        BigNum Tab = 0;
        if (b - a == 1) {
            if (a == 0) {
                Pab = Qab = 1;
            } else {
                Pab = (6 * a - 5) * (2 * a - 1) * (6 * a - 1);
                Qab = a * a * a * C3_OVER_24;
            }
            Tab = Pab * (13591409 + 545140134*a);
            if (a & 1) {
                Tab.setSign(Tab.getSign() * -1);
            }
        } else {
            int64_t m = (a + b) / 2;
            auto res_a = bs(a, m);
            auto res_b = bs(m, b);
            Pab = res_a.first * res_b.first;
            Qab = res_a.second * res_b.second;
            Tab = res_b.second * res_a.third + res_a.first * res_b.third;
        }

        triple res = {Pab, Qab, Tab};
        return res;
    };

    BigNum digits_per_term = "14.181647462725477655525521678181770863769125289828726959816854332";
    int64_t n = static_cast<int>(precision / 14.1816474627 + 1);
    auto res = bs(0, n);

    BigNum one = 1;
    for (int i = 0; i < precision; ++i) {
        one = one * 10;
    }

    auto sqrtC = sqrt(10005 * one, one);

    std::cout << res.second.toString() << '\n';
    std::cout << res.third.toString() << '\n';
//    auto r = res.second * BigNum("1000249968757810059447921878763577780015950243686963146571355115696509678538643042923111879484999732977551938893695661811101310349073901991031130110817620021084773209484713755522300964025814877362597996538778520690844138517442330859859351368356455959879200742054803370888442591728291463927179974654449569355332351563623934371444461273353426110551814049847091601969015554438382830530991787981536628540422125155591789551900389567969478750437233467227464892504600935861206393398280616207431963680961315612113677398957195886501859014838681740375258974396156572830025905031781793481127928508115201563959053331833991286790700407982081441386506614037421393224134385067616132531528652053801465767617231369828928041492395588277347299224032881241981474836093181268398950719719184613180485442334808803907336200724905008814141075979996361648210336419155533604357564303407077666606427135655088765678439754372321736055432944086373962522100785904153733709154299465285532563568099289783119057016422898294516743742803726") * 426880;
    auto r = res.second * sqrtC * 426880;
    std::cout << r.toString() << '\n';
    auto pi = BigNumDiv(r, res.third, 0);
    auto str = pi.toString();
    str.pop_back();
    if (pi.getExp() < 0) {
        for (int i = 0; i < -pi.getExp() + 1; ++i) {
            str.pop_back();
        }
    }
    pi = BigNum(str);
    pi.setExp(-precision + 1);
    return pi;
    return BigNumDiv((1000249968757_bn * 426880_bn), 13591409_bn, 0);
}

// 3.1415926535897932384626433832795028841971693993751058209749445923
// 3.14159265358979323846264338327950288419716939937507832635983263598326359832635983263598326359832636
// 3.1415926535897932384626433832795028841971693993751058209749445923078188
// 3.1415926535897932384626433832795028841971693993751058209749445923078188000