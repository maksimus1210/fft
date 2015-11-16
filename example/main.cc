#include <fft/fft.h>
#include <iostream>

#include <vector>
#include <algorithm>

int main() {
    auto n = 128;
    std::vector<std::complex<double> > values(n);
    auto x = 0.0;
    std::generate(std::begin(values), std::end(values), [&x] { return std::complex<double>(x++, 0.0); });
    fft::fft(std::begin(values), std::end(values));
    fft::fft_inv(std::begin(values), std::end(values));
    for (auto x : values) {
        std::cout << x << std::endl;
    } 
    return 0;    
}
