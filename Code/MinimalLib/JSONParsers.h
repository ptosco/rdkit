#pragma once

#include <GraphMol/FileParsers/PNGParser.h>

namespace RDKit {
namespace MinimalLib {

void updatePropertyPickleOptionsFromJSON(unsigned int &propFlags,
                                         const char *details_json);
void updatePNGMetadataParamsFromJSON(PNGMetadataParams &params,
                                     const char *details_json);

}  // end namespace MinimalLib
}  // end namespace RDKit
