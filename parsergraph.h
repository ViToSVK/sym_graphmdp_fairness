#ifndef _PARSER_H
#define _PARSER_H

#include "streettgraph.h"

class parser {
    
    public:
    
    int basic(streettgraph &G);

    int product(std::string filename, streettgraph &G);
    
    int hoa(std::string filename, streettgraph &G);

};

#endif /* _PARSER_H */
