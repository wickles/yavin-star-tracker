/**
* This work belongs to the Star Tracker project at UCSB.
* Copyright 2010 Alexander Wickes
* Copyright 2010 Gil Tabak
* Copyright 2012 Karanbir Toor
* This work is licensed under the Apache License Version 2.0,
* subject to all terms as reproduced in the included LICENSE file.
*/

#include "ini_reader.h"
#include <cstring>
#include <cstdlib>

static char buffer[256];

int IniGetInt(FILE* file, const char* param, int notfound)
{
    if (!file || !param)
        return notfound;
    rewind(file);
    while(!feof(file))
    {
        fgets(buffer, sizeof(buffer), file);
        if (strncmp(buffer, param, strlen(param)) == 0)
		{
			int r;
			sscanf(buffer, "%*s = %i", &r);
			return r;
        }
    }
    return notfound;
}

float IniGetFloat(FILE* file, const char* param, float notfound)
{
    if (!file || !param)
        return notfound;
    rewind(file);
    while(!feof(file))
    {
        fgets(buffer, sizeof(buffer), file);
        if (strncmp(buffer, param, strlen(param)) == 0)
        {
			float r;
			sscanf(buffer, "%*s = %f", &r);
			return r;
        }
    }
    return notfound;
}

const char* IniGetString(FILE* file, const char* param, char* out, const char* notfound)
{
    if (!out)
        return NULL;
    if (!file || !param)
    {
        if (notfound)
            strcpy(out, notfound);
        return out;
    }
    rewind(file);
    while(!feof(file))
    {
        fgets(buffer, sizeof(buffer), file);
        if (strncmp(buffer, param, strlen(param)) == 0)
        {
            if (sscanf(buffer, "%*s = %s", out) != 1)
            {
                if (notfound)
                    strcpy(out, notfound);
                else
                    out[0] = '\0';
			}
			return out;
        }
	}
	if (notfound)
		strcpy(out, notfound);
	else
		out[0] = '\0';
    return out;
}

bool IniGetBool(FILE* file, const char* param, bool notfound)
{
	if (!file || !param)
		return notfound;
	char buffer[256];
	IniGetString(file, param, buffer, "true");

	if (stricmp(buffer, "true") == 0)
		return true;
	if (stricmp(buffer, "false") == 0)
		return false;
	if (!atoi(buffer))
		return false;

	return true;
}
