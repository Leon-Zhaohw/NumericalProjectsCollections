/*
 * Generate Thue-Morse word prefixes
 * Abdalla Ahmed
 * 2016-01-01
 */


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <assert.h>

typedef std::string                         String;
typedef std::vector<String>                 Strings;

String morph(String &s) {
    String result;
    result.resize(2 * s.length());
    for (int i = 0; i < s.length(); i++) {
        result[2*i] = s[i];
        result[2*i + 1] = (s[i] == '0' ? '1' : '0');
    }
    return result;
}

const char *USAGE_MESSAGE = "Usage: %s <depth>\n";

int main(int argc,char **argv) {
    if (argc != 2) {
        fprintf(stderr, USAGE_MESSAGE, argv[0]);
        exit(1);
    }
    int depth = atoi(argv[1]);
    String s = "0";
    for (int i = 0; i < depth; i++) s = morph(s);
    std::cout << s << "\n";
}