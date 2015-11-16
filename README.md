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
