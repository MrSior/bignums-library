#include <iostream>
#include <string>
#include <concepts>
#include "BigNum.h"

int main()
{
    BigNum a = 11;
    a.printBlocks();
    BigNum b = 1;
    b.printBlocks();
    BigNum c = a * b;
    c.printBlocks();
//    a = b;
    std::cout << (c != a);
    return 0;
}