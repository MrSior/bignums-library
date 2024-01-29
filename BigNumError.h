#ifndef BIGNUMSLIBRARY_BIGNUMERROR_H
#define BIGNUMSLIBRARY_BIGNUMERROR_H

#include <exception>
#include <string>

class BigNumError : public std::exception{
private:
    std::string message_;
public:
    BigNumError(const std::string& ctx) {
        message_ = ctx;
    }

    [[nodiscard]] const char* what() const noexcept override {
        return message_.c_str();
    }
};


#endif //BIGNUMSLIBRARY_BIGNUMERROR_H
