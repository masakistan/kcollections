#pragma once

#include <cstdint>
#include <cstring>

/** 256-bit presence mask (replaces vendored uint256_t). */
struct PrefMask {
  uint64_t w[4];

  PrefMask() { w[0] = w[1] = w[2] = w[3] = 0; }

  explicit PrefMask(uint64_t a, uint64_t b = 0, uint64_t c = 0, uint64_t d = 0) {
    w[0] = a;
    w[1] = b;
    w[2] = c;
    w[3] = d;
  }

  bool test_bit(unsigned bit) const { return (w[bit / 64] >> (bit % 64)) & 1ULL; }

  void set_bit(unsigned bit) { w[bit / 64] |= (1ULL << (bit % 64)); }

  PrefMask operator>>(unsigned shift) const {
    PrefMask out;
    if (shift >= 256) {
      return out;
    }
    unsigned word = shift / 64;
    unsigned off = shift % 64;
    for (int i = 0; i < 4; ++i) {
      if (i + (int)word < 4) {
        out.w[i] = w[i + word] >> off;
        if (off && i + (int)word + 1 < 4) {
          out.w[i] |= w[i + word + 1] << (64 - off);
        }
      }
    }
    return out;
  }

  PrefMask operator<<(unsigned shift) const {
    PrefMask out;
    if (shift >= 256) {
      return out;
    }
    unsigned word = shift / 64;
    unsigned off = shift % 64;
    for (int i = 3; i >= 0; --i) {
      if (i >= (int)word) {
        out.w[i] = w[i - word] << off;
        if (off && i - (int)word - 1 >= 0) {
          out.w[i] |= w[i - word - 1] >> (64 - off);
        }
      }
    }
    return out;
  }

  PrefMask operator|(const PrefMask& other) const {
    PrefMask out;
    for (int i = 0; i < 4; ++i) {
      out.w[i] = w[i] | other.w[i];
    }
    return out;
  }

  PrefMask& operator|=(const PrefMask& other) {
    for (int i = 0; i < 4; ++i) {
      w[i] |= other.w[i];
    }
    return *this;
  }

  static int popcount_vidx(PrefMask vertices, uint8_t bts) {
    vertices = vertices << (256 - (unsigned)bts);
    int vidx = 0;
    for (int i = 0; i < 4; ++i) {
#if defined(__GNUC__) || defined(__clang__)
      vidx += __builtin_popcountll(vertices.w[i]);
#else
      uint64_t x = vertices.w[i];
      while (x) {
        vidx += x & 1;
        x >>= 1;
      }
#endif
    }
    return vidx;
  }
};
