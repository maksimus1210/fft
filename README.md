# fft

## How to include headers?

``` #include <fft/fft.h>  ```

## How to perform Fast Fourier transform?
``` fft::fft(std::begin(values), std::end(values)); ```

## How to perform Inverse Fast Fourier transform?
``` fft::fft(std::begin(values), std::end(values)); ```

## What are the requirements for object stored in variable called `values`?
* std::is_same<decltype(values)::value_type, std::complex<double>>::value == true, 
* std::is_same<std::random_access_iterator_tag,  typename std::iterator_traits<decltype(values)::iterator>::iterator_category>::value == true
* values.size() is power of 2

## What about example?
```
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
```
