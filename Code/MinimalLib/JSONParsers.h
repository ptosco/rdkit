#pragma once

namespace RDKit {
namespace MinimalLib {

void updatePropertyPickleOptionsFromJSON(const char *details_json,
                                         unsigned int &propFlags);

}  // end namespace MinimalLib
}  // end namespace RDKit
