//
//  Tokenizer.hpp
//

#ifndef Tokenizer_hpp
#define Tokenizer_hpp

#include <string>
#include <vector>

/* -------------------------
 
 CLASS TOKENIZER
 
 --------------------------- */

const std::string DEFAULT_DELIMITER = " \t\v\n\r\f";

class Tokenizer
{
public:
    // ctor/dtor
    Tokenizer();
    Tokenizer(const std::string& str, const std::string& delimiter=DEFAULT_DELIMITER);
    ~Tokenizer();
    
    // set string and delimiter
    void set(const std::string& str, const std::string& delimiter=DEFAULT_DELIMITER);
    void setString(const std::string& str);             // set source string only
    void setDelimiter(const std::string& delimiter);    // set delimiter string only
    
    std::string next();                                 // return the next token, return "" if it ends
    
    std::vector<std::string> split();                   // return array of tokens from current cursor
    
protected:
    
private:
    void skipDelimiter();                               // ignore leading delimiters
    bool isDelimiter(char c);                           // check if the current char is delimiter
    
    std::string buffer;                                 // input string
    std::string token;                                  // output string
    std::string delimiter;                              // delimiter string
    std::string::const_iterator currPos;                // string iterator pointing the current position
    
};


#endif 
