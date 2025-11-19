#pragma once
#include <cmath>
#include <limits>
#include <type_traits>
#include <ostream>

template <class T>
struct WelfordStats {
    static_assert(std::is_arithmetic<T>::value, "T must be arithmetic");
    static_assert(std::is_standard_layout<T>::value, "T must be standard layout");
    static_assert(std::is_trivially_copyable<T>::value, "T must be trivially copyable");
    static_assert(std::is_trivial<T>::value, "T must be trivial");

    // POD fields — trivially copyable
    long long n;
    long double mean;
    long double M2;
    T min_val;
    T max_val;

    // Default constructor
    constexpr WelfordStats()
        : n(0),
          mean(0),
          M2(0),
          min_val(std::numeric_limits<T>::max()),
          max_val(std::numeric_limits<T>::lowest())
    {}

    // Merging constructor: merges b into a
    constexpr WelfordStats(const WelfordStats& a, const WelfordStats& b)
        : WelfordStats(a) {
        this->merge(b);
    }

    // Reset to initial state
    void reset() {
        n = 0;
        mean = 0;
        M2 = 0;
        min_val = std::numeric_limits<T>::max();
        max_val = std::numeric_limits<T>::lowest();
    }

    // Add one value
    void add(T x) {
        if (x < min_val) min_val = x;
        if (x > max_val) max_val = x;

        n++;
        long double delta = x - mean;
        mean += delta / n;
        long double delta2 = x - mean;
        M2 += delta * delta2;
    }

    // Merge another stats object
    void merge(const WelfordStats& other) {
        if (other.n == 0) return;
        if (n == 0) { *this = other; return; }

        if (other.min_val < min_val) min_val = other.min_val;
        if (other.max_val > max_val) max_val = other.max_val;

        long double n_total = n + other.n;
        long double delta = other.mean - mean;

        mean += delta * (other.n / n_total);

        M2 += other.M2 +
              delta * delta * (n * other.n / n_total);

        n = static_cast<long long>(n_total);
    }

    // += T
    WelfordStats &operator+=(T x) {
        add(x);
        return *this;
    }

    // += WelfordStats
    WelfordStats &operator+=(const WelfordStats& other) {
        merge(other);
        return *this;
    }

    // + T (non-mutating)
    WelfordStats operator+(T x) const {
        WelfordStats tmp = *this;
        tmp.add(x);
        return tmp;
    }

    // + WelfordStats (non-mutating merge)
    WelfordStats operator+(const WelfordStats& other) const {
        return WelfordStats(*this, other);
    }

    // *= long double
    WelfordStats &operator*=(long double k) {
        mean *= k;
        M2   *= (k * k);
        if (n > 0) {
            min_val = static_cast<T>(min_val * k);
            max_val = static_cast<T>(max_val * k);
        }
        return *this;
    }

    // * long double (non-mutating scale)
    WelfordStats operator*(long double k) const {
        WelfordStats tmp = *this;
        tmp *= k;
        return tmp;
    }

    // Accessors
    bool empty() const { return n == 0; }
    T min() const { return min_val; }
    T max() const { return max_val; }
    long double count() const { return n; }
    long double mean_value() const { return mean; }

    long double variance_population() const {
        return (n > 0) ? M2 / n : 0;
    }

    long double variance_sample() const {
        return (n > 1) ? M2 / (n - 1) : 0;
    }

    long double stddev() const {
        return std::sqrt(variance_population());
    }
};

// Print helper
template <class T>
std::ostream &operator<<(std::ostream& os, const WelfordStats<T> &s) {
    os << "WelfordStats{"
       << "n=" << s.count()
       << ", min=" << s.min()
       << ", max=" << s.max()
       << ", mean=" << s.mean_value()
       << ", variance=" << s.variance_population()
       << ", stddev=" << s.stddev()
       << "}";
    return os;
}
