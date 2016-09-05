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

// Comment out to disable SDL
//#define SDL_ENABLED

#include <windows.h>
#include "BUF_USBCCDCamera_SDK.h"
#ifdef SDL_ENABLED
#include "SDL/SDL.h"
#include "SDL/SDL_ttf.h"
#endif

using namespace std;

#endif // !defined(AFX_STDAFX_H__B6EAA554_08E3_488D_A3D5_DA3AE5370FA8__INCLUDED_)