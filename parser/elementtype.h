#ifndef __ELEMENTTYPE_H
#define __ELEMENTTYPE_H

// types of elements

enum class ElementType {
  _NONE = -1,
  _MARKER = 1,
  _DRIFT = 2,
  _RF = 3,
  _SBEND = 4, 
  _QUAD  = 5,
  _SEXTUPOLE  = 6,
  _OCTUPOLE = 7,
  _MULT  = 8,
  _SOLENOID = 9,
  _ELEMENT = 10,
  _LINE = 11,
  _REV_LINE= -11, //for line inversion in sublines
  _COLLIMATOR = 12, // obsolete?
  _ECOL = 13,
  _MUSPOILER = 62,
  _RCOL = 14,
  _LASER=15,
  _MATERIAL=16,
  _RBEND=17,
  _ATOM = 18,
  _SEQUENCE = 19,
  _SCREEN = 21,
  _AWAKESCREEN = 22,
  _VKICK=31,
  _HKICK=32,
  _SAMPLER = 41,
  _CSAMPLER = 42,
  _DUMP = 43,
  _GAS = 51,
  _TUNNEL = 52,
  _TRANSFORM3D = 61,
  _TELEPORTER  = 98,
  _TERMINATOR  = 99
};

namespace GMAD {
  const char *typestr(ElementType type);
}

#endif
