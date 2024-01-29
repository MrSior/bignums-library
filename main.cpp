#include <iostream>
#include <string>
#include <concepts>
#include "BigNum.h"

int main()
{
    BigNum a = 11;
    a.printBlocks();
    BigNum b = 1000;
    b.printBlocks();
    BigNum c = a * b;
    c.printBlocks();
//    a = b;
    std::cout << c.toString();
    return 0;
}