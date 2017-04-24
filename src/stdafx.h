// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#if !defined(AFX_STDAFX_H__B6EAA554_08E3_488D_A3D5_DA3AE5370FA8__INCLUDED_)
#define AFX_STDAFX_H__B6EAA554_08E3_488D_A3D5_DA3AE5370FA8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

// TODO: reference additional headers your program requires here

// Moved to project configuration
//#define SDL_ENABLED

#include <windows.h>
#include "BUF_USBCCDCamera_SDK.h"
#ifdef SDL_ENABLED
#include "SDL/SDL.h"
#include "SDL/SDL_ttf.h"
#endif

// fix for VS-2015+
// also need to define _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS
// to allow use of deprecated container hash_map
#if (_MSC_VER >= 1900) && defined(SDL_ENABLED)
FILE _iob[] = { *stdin, *stdout, *stderr };
extern "C" FILE * __cdecl __iob_func(void)
{
	return _iob;
}
#endif

//using namespace std;

#endif // !defined(AFX_STDAFX_H__B6EAA554_08E3_488D_A3D5_DA3AE5370FA8__INCLUDED_)