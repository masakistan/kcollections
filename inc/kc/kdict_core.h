#pragma once

#ifndef KC_KDICT_CORE_INCLUDED
#define KC_KDICT_CORE_INCLUDED

#define KC_CONTAINER KcDictContainer
#define KC_VERTEX KcDictVertex
#define KC_UC KcDictUC
#define KC_ITERATOR KcDictIterator

#define KDICT 1

#include "globals.h"
#include "helper.h"
#include "kc_io.h"
#include "pref_mask.h"
#include "UContainer.h"
#include "Vertex.h"
#include "Kcontainer.h"

namespace kc::dict {
template <typename T>
using Container = KcDictContainer<T>;
template <typename T>
using Vertex = KcDictVertex<T>;
template <typename T>
using UC = KcDictUC<T>;
}  // namespace kc::dict

#endif
