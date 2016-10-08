//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_OPENMM_H__
#define __RD_OPENMM_H__

#include <OpenMM.h>

namespace OpenMM {
  class OpenMMPlugins {
    public:
      static void loadOpenMMPlugins(
          const boost::uint8_t *aromatic_types = NULL);
      
}

#endif
