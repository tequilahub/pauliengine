#include "PauliString.h"

int main() {
    std::complex<double> i{0.0, 1.0};
    PauliString<std::complex<double>> w0(std::vector<uint64_t>{0}, std::vector<uint64_t>{0}, 1);
    PauliString<std::complex<double>> x0(std::unordered_map<int, std::string>{{0, "X"}}, 1);
    PauliString<std::complex<double>> y0(std::unordered_map<int, std::string>{{0, "Y"}}, 1);
    PauliString<std::complex<double>> z0(std::unordered_map<int, std::string>{{0, "Z"}}, 1);
    bool pass = true;

    // Check single Pauli multiplication table.
    pass &= w0 * w0 == w0;
    pass &= w0 * x0 == x0;
    pass &= w0 * y0 == y0;
    pass &= w0 * z0 == z0;
    pass &= x0 * w0 == x0;
    pass &= x0 * x0 == w0;
    pass &= x0 * y0 == z0 * i;
    pass &= x0 * z0 == y0 * -i;
    pass &= y0 * w0 == y0;
    pass &= y0 * x0 == z0 * -i;
    pass &= y0 * y0 == w0;
    pass &= y0 * z0 == x0 * i;
    pass &= z0 * w0 == z0;
    pass &= z0 * x0 == y0 * -i;
    pass &= z0 * y0 == x0 * i;
    pass &= z0 * z0 == w0;

    // Check automatic length extension.
    PauliString<std::complex<double>> v1(std::unordered_map<int, std::string>{{50, "X"}}, 1);
    PauliString<std::complex<double>> v2(std::unordered_map<int, std::string>{{100, "Z"}}, 1);
    PauliString<std::complex<double>> v12(std::unordered_map<int, std::string>{{50, "X"}, {100, "Z"}}, 1);
    pass &= v1.x.size() == 1;
    pass &= v2.x.size() == 2;
    pass &= (v1 * v2).x.size() == 2;
    pass &= v1 * v2 == v12;
    pass &= v2 * v1 == v12;
    v1 *= v2;
    pass &= v1 == v12;
    pass &= v1.x.size() == 2;

    // Check accumulation mod 4 within a word.
    {
        PauliString<std::complex<double>> xxxxx3(std::unordered_map<int, std::string>{{0, "X"}, {1, "X"}, {2, "X"}, {3, "X"}, {4, "X"}}, 3);
        PauliString<std::complex<double>> z5(std::unordered_map<int, std::string>{{0, "Z"}}, 5);
        PauliString<std::complex<double>> zz5(std::unordered_map<int, std::string>{{0, "Z"}, {1, "Z"}}, 5);
        PauliString<std::complex<double>> zzz5(std::unordered_map<int, std::string>{{0, "Z"}, {1, "Z"}, {2, "Z"}}, 5);
        PauliString<std::complex<double>> zzzz5(std::unordered_map<int, std::string>{{0, "Z"}, {1, "Z"}, {2, "Z"}, {3, "Z"}}, 5);
        PauliString<std::complex<double>> zzzzz5(std::unordered_map<int, std::string>{{0, "Z"}, {1, "Z"}, {2, "Z"}, {3, "Z"}, {4, "Z"}}, 5);
        pass &= xxxxx3 * z5 == PauliString<std::complex<double>>(std::unordered_map<int, std::string>{{0, "Y"}, {1, "X"}, {2, "X"}, {3, "X"}, {4, "X"}}, 15.0 * -i);
        pass &= xxxxx3 * zz5 == PauliString<std::complex<double>>(std::unordered_map<int, std::string>{{0, "Y"}, {1, "Y"}, {2, "X"}, {3, "X"}, {4, "X"}}, 15.0 * -1.0);
        pass &= xxxxx3 * zzz5 == PauliString<std::complex<double>>(std::unordered_map<int, std::string>{{0, "Y"}, {1, "Y"}, {2, "Y"}, {3, "X"}, {4, "X"}}, 15.0 * i);
        pass &= xxxxx3 * zzzz5 == PauliString<std::complex<double>>(std::unordered_map<int, std::string>{{0, "Y"}, {1, "Y"}, {2, "Y"}, {3, "Y"}, {4, "X"}}, 15.0);
        pass &= xxxxx3 * zzzzz5 == PauliString<std::complex<double>>(std::unordered_map<int, std::string>{{0, "Y"}, {1, "Y"}, {2, "Y"}, {3, "Y"}, {4, "Y"}}, 15.0 * -i);
    }

    // Check accumulation mod 4 across words.
    {
        PauliString<std::complex<double>> xxxxx3(std::unordered_map<int, std::string>{{0, "X"}, {100, "X"}, {200, "X"}, {300, "X"}, {400, "X"}}, 3);
        PauliString<std::complex<double>> z5(std::unordered_map<int, std::string>{{0, "Z"}}, 5);
        PauliString<std::complex<double>> zz5(std::unordered_map<int, std::string>{{0, "Z"}, {100, "Z"}}, 5);
        PauliString<std::complex<double>> zzz5(std::unordered_map<int, std::string>{{0, "Z"}, {100, "Z"}, {200, "Z"}}, 5);
        PauliString<std::complex<double>> zzzz5(std::unordered_map<int, std::string>{{0, "Z"}, {100, "Z"}, {200, "Z"}, {300, "Z"}}, 5);
        PauliString<std::complex<double>> zzzzz5(std::unordered_map<int, std::string>{{0, "Z"}, {100, "Z"}, {200, "Z"}, {300, "Z"}, {400, "Z"}}, 5);
        pass &= xxxxx3 * z5 == PauliString<std::complex<double>>(std::unordered_map<int, std::string>{{0, "Y"}, {100, "X"}, {200, "X"}, {300, "X"}, {400, "X"}}, 15.0 * -i);
        pass &= xxxxx3 * zz5 == PauliString<std::complex<double>>(std::unordered_map<int, std::string>{{0, "Y"}, {100, "Y"}, {200, "X"}, {300, "X"}, {400, "X"}}, 15.0 * -1.0);
        pass &= xxxxx3 * zzz5 == PauliString<std::complex<double>>(std::unordered_map<int, std::string>{{0, "Y"}, {100, "Y"}, {200, "Y"}, {300, "X"}, {400, "X"}}, 15.0 * i);
        pass &= xxxxx3 * zzzz5 == PauliString<std::complex<double>>(std::unordered_map<int, std::string>{{0, "Y"}, {100, "Y"}, {200, "Y"}, {300, "Y"}, {400, "X"}}, 15.0);
        pass &= xxxxx3 * zzzzz5 == PauliString<std::complex<double>>(std::unordered_map<int, std::string>{{0, "Y"}, {100, "Y"}, {200, "Y"}, {300, "Y"}, {400, "Y"}}, 15.0 * -i);
    }

    if (pass) {
        std::cerr << "passed\n";
    } else {
        std::cerr << "FAILED\n";
    }
    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
