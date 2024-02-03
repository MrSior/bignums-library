#include <iostream>
#include <string>
#include <concepts>
#include "BigNum.h"

using namespace std::chrono;

int main()
{
//    BigNum a = "000000166350886035.37698900000";
////    BigNum a = "0.192347";
//    a.printBlocks();
//    BigNum b = "864847832487";
//    b.printBlocks();
//
//    std::cout << (a / b).toString();

//    BigNum a = 1_bn;
//    a.printBlocks();
//    BigNum b = 239;
//    b.printBlocks();
//
//
//    std::cout << BigNumDiv(a, b, 100, false).toString() << std::endl;
//
//    b = b * b;
//
//    BigNum c = CalcPi(10);
//    c.printBlocks();
//
//    BigNum d = 10_bn;
//
////    Division(a, b, 1);
//
////    a = "7986610278667.6989";
////    b = "864847832487";
////    BigNum c = a - 9 * b;
////    c.printBlocks();
//
//    std::cout << "res " << c.toString() << std::endl;
////    std::cout << (c <= d) << std::endl;
//

//    auto arctan = [](const BigNum& x, int steps) {
//        BigNum res = x;
//        auto numer = x * x * x;
//        int denom = 3;
//        std::cout << x.toString() << std::endl;
//        for (int i = 1; i < steps; ++i, denom += 2, numer = numer * x * x) {
//            auto sum = BigNumDiv(numer, denom, 20);
//            std::cout << numer.toString() << " / " << denom << " = " << sum.toString() << std::endl;
//            if (i % 2 == 0) {
//                res = res + sum;
//            } else {
//                res = res - sum;
//            }
//        }
//        return res;
//    };

    auto time_s = high_resolution_clock::now();
    std::cout << CalcPi(0).toString() << std::endl;
    auto time_end = high_resolution_clock::now();
    std::cout << "Time: " << duration_cast<seconds>(time_end - time_s).count() << " seconds" << std::endl;

    return 0;
}