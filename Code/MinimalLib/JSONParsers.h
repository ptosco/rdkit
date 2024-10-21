#pragma once

namespace RDKit {
namespace MinimalLib {

bool updatePropertyPickleOptionsFromJSON(const char *details_json,
                                         unsigned int &propFlags);

}  // end namespace MinimalLib
}  // end namespace RDKit
