#ifndef TRUSS_FEM_JSON_PARSER_H
#define TRUSS_FEM_JSON_PARSER_H

#include <string>
#include "TrussModel.h"

namespace truss_fem {

class JsonParser {
public:
    static TrussModel parse(const std::string& filepath);
};

} // namespace truss_fem

#endif
