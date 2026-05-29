#pragma once

#ifndef KC_KCOUNTER_CORE_INCLUDED
#define KC_KCOUNTER_CORE_INCLUDED

#define KC_CONTAINER KcCounterContainer
#define KC_VERTEX KcCounterVertex
#define KC_UC KcCounterUC
#define KC_ITERATOR KcCounterIterator

#define KCOUNTER 1

#include "globals.h"
#include "helper.h"
#include "kc_io.h"
#include "pref_mask.h"
#include "UContainer.h"
#include "Vertex.h"
#include "Kcontainer.h"

namespace kc::counter {
using Container = KcCounterContainer<int>;
using Vertex = KcCounterVertex<int>;
using UC = KcCounterUC<int>;
}  // namespace kc::counter

#endif
