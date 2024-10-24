#define USE_BETTER_ENUMS
#include "JSONParsers.h"
#include <GraphMol/MolPickler.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace MinimalLib {

void updatePropertyPickleOptionsFromJSON(unsigned int &propFlags,
                                         const char *details_json) {
  if (details_json && strlen(details_json)) {
    std::istringstream ss;
    boost::property_tree::ptree pt;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    const auto nodeIt = pt.find("propertyFlags");
    if (nodeIt != pt.not_found()) {
      auto propertyFlagsFromJson =
          (+PicklerOps::PropertyPickleOptions::NoProps)._to_integral();
      for (const auto *key : PicklerOps::PropertyPickleOptions::_names()) {
        if (nodeIt->second.get(key, false)) {
          propertyFlagsFromJson |=
              PicklerOps::PropertyPickleOptions::_from_string(key)
                  ._to_integral();
        }
      }
      propFlags = propertyFlagsFromJson;
    }
  }
}

}  // end namespace MinimalLib
}  // end namespace RDKit
