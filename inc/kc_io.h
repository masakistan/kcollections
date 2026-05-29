#pragma once

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "pref_mask.h"

namespace kc_io {

constexpr char FILE_MAGIC[4] = {'K', 'C', 'O', 'L'};
constexpr uint16_t FILE_VERSION = 2;

enum class ContainerKind : uint8_t {
  KIND_SET = 1,
  KIND_DICT = 2,
  KIND_COUNTER = 3,
};

namespace detail {

inline void write_u8(std::ostream& os, uint8_t v) { os.put(static_cast<char>(v)); }

inline void write_u16_le(std::ostream& os, uint16_t v) {
  write_u8(os, static_cast<uint8_t>(v & 0xff));
  write_u8(os, static_cast<uint8_t>((v >> 8) & 0xff));
}

inline void write_u32_le(std::ostream& os, uint32_t v) {
  for (int i = 0; i < 4; ++i) {
    write_u8(os, static_cast<uint8_t>((v >> (8 * i)) & 0xff));
  }
}

inline void write_u64_le(std::ostream& os, uint64_t v) {
  for (int i = 0; i < 8; ++i) {
    write_u8(os, static_cast<uint8_t>((v >> (8 * i)) & 0xff));
  }
}

inline uint8_t read_u8(std::istream& is) {
  int c = is.get();
  if (c == EOF) {
    throw std::runtime_error("kcollections: failed to read archive (truncated)");
  }
  return static_cast<uint8_t>(c);
}

inline uint16_t read_u16_le(std::istream& is) {
  uint16_t v = read_u8(is);
  v |= static_cast<uint16_t>(read_u8(is)) << 8;
  return v;
}

inline uint32_t read_u32_le(std::istream& is) {
  uint32_t v = 0;
  for (int i = 0; i < 4; ++i) {
    v |= static_cast<uint32_t>(read_u8(is)) << (8 * i);
  }
  return v;
}

inline uint64_t read_u64_le(std::istream& is) {
  uint64_t v = 0;
  for (int i = 0; i < 8; ++i) {
    v |= static_cast<uint64_t>(read_u8(is)) << (8 * i);
  }
  return v;
}

template <typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type
write_arithmetic_le(std::ostream& os, T value) {
  if (sizeof(T) == sizeof(uint32_t)) {
    uint32_t bits = 0;
    static_assert(sizeof(T) == sizeof(bits), "unexpected float size");
    std::memcpy(&bits, &value, sizeof(T));
    write_u32_le(os, bits);
  } else {
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(T));
    write_u64_le(os, bits);
  }
}

template <typename T>
inline typename std::enable_if<!std::is_floating_point<T>::value, void>::type
write_arithmetic_le(std::ostream& os, T value) {
  static_assert(std::is_arithmetic<T>::value, "expected arithmetic type");
  typename std::make_unsigned<typename std::decay<T>::type>::type bits =
      static_cast<typename std::make_unsigned<typename std::decay<T>::type>::type>(value);
  for (size_t i = 0; i < sizeof(T); ++i) {
    write_u8(os, static_cast<uint8_t>((bits >> (8 * i)) & 0xff));
  }
}

template <typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type
read_arithmetic_le(std::istream& is, T& value) {
  if (sizeof(T) == sizeof(uint32_t)) {
    uint32_t bits = read_u32_le(is);
    std::memcpy(&value, &bits, sizeof(T));
  } else {
    uint64_t bits = read_u64_le(is);
    std::memcpy(&value, &bits, sizeof(T));
  }
}

template <typename T>
inline typename std::enable_if<!std::is_floating_point<T>::value, void>::type
read_arithmetic_le(std::istream& is, T& value) {
  static_assert(std::is_arithmetic<T>::value, "expected arithmetic type");
  typename std::make_unsigned<typename std::decay<T>::type>::type bits = 0;
  for (size_t i = 0; i < sizeof(T); ++i) {
    bits |= static_cast<typename std::make_unsigned<typename std::decay<T>::type>::type>(read_u8(is))
            << (8 * i);
  }
  value = static_cast<T>(bits);
}

}  // namespace detail

class BinaryOutArchive {
  std::ostream& os_;

public:
  explicit BinaryOutArchive(std::ostream& os) : os_(os) {}

  template <typename T>
  typename std::enable_if<std::is_arithmetic<typename std::decay<T>::type>::value,
                          BinaryOutArchive&>::type
  operator&(const T& value) {
    detail::write_arithmetic_le(os_, value);
    check();
    return *this;
  }

  BinaryOutArchive& operator&(const std::string& value) {
    const uint64_t len = value.size();
    (*this) & len;
    if (len > 0) {
      os_.write(value.data(), static_cast<std::streamsize>(len));
      check();
    }
    return *this;
  }

  BinaryOutArchive& operator&(const PrefMask& value) {
    for (int i = 0; i < 4; ++i) {
      (*this) & value.w[i];
    }
    return *this;
  }

  template <typename T>
  BinaryOutArchive& operator&(const std::vector<T>& value) {
    const uint64_t len = value.size();
    (*this) & len;
    for (const T& item : value) {
      (*this) & item;
    }
    return *this;
  }

private:
  void check() {
    if (!os_) {
      throw std::runtime_error("kcollections: failed to write archive");
    }
  }
};

class BinaryInArchive {
  std::istream& is_;

public:
  explicit BinaryInArchive(std::istream& is) : is_(is) {}

  template <typename T>
  typename std::enable_if<std::is_arithmetic<typename std::decay<T>::type>::value,
                          BinaryInArchive&>::type
  operator&(T& value) {
    detail::read_arithmetic_le(is_, value);
    return *this;
  }

  BinaryInArchive& operator&(std::string& value) {
    uint64_t len = 0;
    (*this) & len;
    value.resize(static_cast<size_t>(len));
    if (len > 0) {
      is_.read(&value[0], static_cast<std::streamsize>(len));
      if (!is_) {
        throw std::runtime_error("kcollections: failed to read archive (truncated)");
      }
    }
    return *this;
  }

  BinaryInArchive& operator&(PrefMask& value) {
    for (int i = 0; i < 4; ++i) {
      (*this) & value.w[i];
    }
    return *this;
  }

  template <typename T>
  BinaryInArchive& operator&(std::vector<T>& value) {
    uint64_t len = 0;
    (*this) & len;
    value.resize(static_cast<size_t>(len));
    for (T& item : value) {
      (*this) & item;
    }
    return *this;
  }
};

inline void write_file_header(BinaryOutArchive& ar, ContainerKind kind) {
  ar & FILE_MAGIC[0] & FILE_MAGIC[1] & FILE_MAGIC[2] & FILE_MAGIC[3];
  const uint16_t version = FILE_VERSION;
  const uint8_t kind_byte = static_cast<uint8_t>(kind);
  ar & version;
  ar & kind_byte;
}

inline void read_file_header(BinaryInArchive& ar, ContainerKind expected) {
  char magic[4] = {};
  ar & magic[0] & magic[1] & magic[2] & magic[3];
  if (magic[0] != FILE_MAGIC[0] || magic[1] != FILE_MAGIC[1] ||
      magic[2] != FILE_MAGIC[2] || magic[3] != FILE_MAGIC[3]) {
    throw std::runtime_error("kcollections: invalid file magic (not a kcollections archive)");
  }
  uint16_t version = 0;
  uint8_t kind_byte = 0;
  ar & version;
  ar & kind_byte;
  if (version != FILE_VERSION) {
    throw std::runtime_error(
        "kcollections: unsupported archive version (expected v2, rebuild index)");
  }
  if (static_cast<ContainerKind>(kind_byte) != expected) {
    throw std::runtime_error("kcollections: archive container type mismatch");
  }
}

}  // namespace kc_io
