/*
 * Utilities for handling words
 * Abdalla Ahmed
 * 2015-12-29
 */

#ifndef WORDS_H
#define WORDS_H


#include <string>
#include <algorithm>

typedef std::string                         String;
typedef std::vector<String>                 Strings;

String neg(const String &S) {                                                   // Replace 0's by 1's and 1's by 0's.
    String result;
    result.resize(S.length());
    for (int i = 0; i < S.length(); i++)
        result[i] = (S[i] == '0' ? '1' : '0');
    return result;
}

String mirror(const String &S) {                                                // Mirror a string left-to-right
    String result;
    int n = S.length();
    result.resize(n);
    for (int i = 0; i < n; i++)
        result[(n-1) - i] = S[i];
    return result;
}

int findOrInsert(const String &s, Strings &list) {
    int i;
    for (i = 0; i < list.size(); i++) {
        if (s.compare( list[i] ) == 0) return i;
    }
    list.push_back(s);
    return i;
}

int findOrAbort(const String &s, const Strings &list) {
    int i;
    for (i = 0; i < list.size(); i++) {
        if (s.compare( list[i] ) == 0) return i;
    }
    fprintf(stderr, "String %s not fount in supplied list\n");
    exit(1);
}

Strings factorize(String s, int width) {
    Strings list;
    for (int i = 0; i < s.length() - width; i++) {
        String sub = s.substr(i, width);
        findOrInsert(sub, list);
    }
    return list;
}

String THMorph(const String &s) {                                               // Apply Thue-Morse morphism to string s
    String result;
    result.resize(2 * s.length());
    for (int i = 0; i < s.length(); i++) {
        result[2*i] = s[i];
        result[2*i + 1] = (s[i] == '0' ? '1' : '0');
    }
    return result;
}

String THMorph(const String &s, int times) {
    String result = s;
    for (int i = 0; i < times; i++) {
        result = THMorph(result);
    }
    return result;
}


#endif
