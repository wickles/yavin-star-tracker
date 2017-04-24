/**
* This work belongs to the Star Tracker project at UCSB.
* Copyright 2010 Alexander Wickes
* Copyright 2010 Gil Tabak
* Copyright 2012 Karanbir Toor
* This work is licensed under the Apache License Version 2.0,
* subject to all terms as reproduced in the included LICENSE file.
*/

#include <cstdio>

int IniGetInt(FILE* file, const char* param, int notfound);

float IniGetFloat(FILE* file, const char* param, float notfound);

const char* IniGetString(FILE* file, const char* param, char* out, const char* notfound);

bool IniGetBool(FILE* file, const char* param, bool notfound);