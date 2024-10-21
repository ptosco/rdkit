#pragma once

#include <GraphMol/FileParsers/PNGParser.h>

namespace RDKit {
namespace MinimalLib {

bool updatePropertyPickleOptionsFromJSON(const char *details_json,
                                         unsigned int &propFlags);
void updatePNGMetadataParamsFromJSON(const char *details_json,
                                     PNGMetadataParams &params);

}  // end namespace MinimalLib
}  // end namespace RDKit
