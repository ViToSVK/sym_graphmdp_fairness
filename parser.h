#ifndef PARSER_H
#define PARSER_H

#include "mec.h"

class parser {
    
public:
    
    int basic(mec &M);

    int product(std::string filename, mec &M, int Rpercent);
    
};

#endif /* PARSER_H */

