#pragma once

#ifndef KC_KSET_CORE_INCLUDED
#define KC_KSET_CORE_INCLUDED

#define KC_CONTAINER KcSetContainer
#define KC_VERTEX KcSetVertex
#define KC_UC KcSetUC
#define KC_ITERATOR KcSetIterator

#define KSET 1

#include "globals.h"
#include "helper.h"
#include "kc_io.h"
#include "pref_mask.h"
#include "UContainer.h"
#include "Vertex.h"
#include "Kcontainer.h"

namespace kc::set {
using Container = KcSetContainer;
using Vertex = KcSetVertex;
using UC = KcSetUC;
}  // namespace kc::set

#endif
