#define USE_BETTER_ENUMS
#include "JSONParsers.h"
#include <GraphMol/MolPickler.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace MinimalLib {

bool updatePropertyPickleOptionsFromJSON(const char *details_json,
                                         unsigned int &propFlags) {
  auto propertyPickleOptionsFromJson =
      (+PicklerOps::PropertyPickleOptions::NoProps)._to_integral();
  bool res = false;
  if (details_json && strlen(details_json)) {
    std::istringstream ss;
    boost::property_tree::ptree pt;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    const auto nodeIt = pt.find("propertyPickleOptions");
    if (nodeIt != pt.not_found()) {
      for (const auto *key : PicklerOps::PropertyPickleOptions::_names()) {
        propertyPickleOptionsFromJson |=
            (nodeIt->second.get(key, false)
                 ? PicklerOps::PropertyPickleOptions::_from_string(key)
                 : +PicklerOps::PropertyPickleOptions::NoProps)
                ._to_integral();
      }
      propFlags = propertyPickleOptionsFromJson;
      res = true;
    }
  }
  return res;
}

}  // end namespace MinimalLib
}  // end namespace RDKit
