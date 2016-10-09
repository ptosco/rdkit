//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//  Copyright (C) 2016 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#if 0
#ifndef __RD_OPENMMHELPERS_H__
#define __RD_OPENMMHELPERS_H__

#include <string>
#include <vector>

namespace RDKit {
namespace OpenMM {

class OpenMMPlugins {
  public:
    const std::vector<std::string>& load(
      const std::string &dir = std::string());
    static OpenMMPlugins *instance();
  private:
    //! to force this to be a singleton, the constructor must be private
    OpenMMPlugins() :
      d_alreadyLoaded(false) {};
    static class OpenMMPlugins *ds_instance;  //!< the singleton
    bool d_alreadyLoaded;
    std::vector<std::string> d_loadedPlugins;
};

} // end namespace OpenMM
} // end namespace RDKit

#endif // __RD_OPENMMHELPERS_H__
#endif // RDK_BUILD_WITH_OPENMM
