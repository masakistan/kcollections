#pragma once

#include <cstdint>
#include <fstream>
#include <list>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "uint128_t.h"
#include "uint256_t.h"

namespace kc_io {

constexpr char FILE_MAGIC[4] = {'K', 'C', 'O', 'L'};
constexpr uint16_t FILE_VERSION = 1;

enum class ContainerKind : uint8_t {
  KIND_SET = 1,
  KIND_DICT = 2,
  KIND_COUNTER = 3,
};

class BinaryOutArchive {
  std::ostream& os_;

public:
  explicit BinaryOutArchive(std::ostream& os) : os_(os) {}

  template <typename T>
  typename std::enable_if<std::is_arithmetic<typename std::decay<T>::type>::value,
                          BinaryOutArchive&>::type
  operator&(const T& value) {
    write_raw(&value, sizeof(T));
    return *this;
  }

  BinaryOutArchive& operator&(const std::string& value) {
    const uint64_t len = value.size();
    (*this) & len;
    if (len > 0) {
      write_raw(value.data(), len);
    }
    return *this;
  }

  BinaryOutArchive& operator&(const uint128_t& value) {
    (*this) & value.UPPER;
    (*this) & value.LOWER;
    return *this;
  }

  BinaryOutArchive& operator&(const uint256_t& value) {
    (*this) & value.UPPER;
    (*this) & value.LOWER;
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

  template <typename T>
  BinaryOutArchive& operator&(const std::set<T>& value) {
    const uint64_t len = value.size();
    (*this) & len;
    for (const T& item : value) {
      (*this) & item;
    }
    return *this;
  }

  template <typename T>
  BinaryOutArchive& operator&(const std::list<T>& value) {
    const uint64_t len = value.size();
    (*this) & len;
    for (const T& item : value) {
      (*this) & item;
    }
    return *this;
  }

private:
  void write_raw(const void* data, size_t nbytes) {
    os_.write(static_cast<const char*>(data), static_cast<std::streamsize>(nbytes));
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
    read_raw(&value, sizeof(T));
    return *this;
  }

  BinaryInArchive& operator&(std::string& value) {
    uint64_t len = 0;
    (*this) & len;
    value.resize(static_cast<size_t>(len));
    if (len > 0) {
      read_raw(&value[0], len);
    }
    return *this;
  }

  BinaryInArchive& operator&(uint128_t& value) {
    (*this) & value.UPPER;
    (*this) & value.LOWER;
    return *this;
  }

  BinaryInArchive& operator&(uint256_t& value) {
    (*this) & value.UPPER;
    (*this) & value.LOWER;
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

  template <typename T>
  BinaryInArchive& operator&(std::set<T>& value) {
    uint64_t len = 0;
    (*this) & len;
    value.clear();
    for (uint64_t i = 0; i < len; ++i) {
      T item{};
      (*this) & item;
      value.insert(std::move(item));
    }
    return *this;
  }

  template <typename T>
  BinaryInArchive& operator&(std::list<T>& value) {
    uint64_t len = 0;
    (*this) & len;
    value.clear();
    for (uint64_t i = 0; i < len; ++i) {
      T item{};
      (*this) & item;
      value.push_back(std::move(item));
    }
    return *this;
  }

private:
  void read_raw(void* data, size_t nbytes) {
    is_.read(static_cast<char*>(data), static_cast<std::streamsize>(nbytes));
    if (!is_) {
      throw std::runtime_error("kcollections: failed to read archive (truncated or invalid file)");
    }
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
    throw std::runtime_error("kcollections: unsupported archive version");
  }
  if (static_cast<ContainerKind>(kind_byte) != expected) {
    throw std::runtime_error("kcollections: archive container type mismatch");
  }
}

}  // namespace kc_io
