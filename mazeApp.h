//---------------------------------------------------------------------------
//
// Name:        mazeApp.h
// Author:      Administrator
// Created:     2/5/2008 9:20:51 PM
// Description: 
//
//---------------------------------------------------------------------------

#ifndef __MAZEDLGApp_h__
#define __MAZEDLGApp_h__

#ifdef __BORLANDC__
	#pragma hdrstop
#endif

#ifndef WX_PRECOMP
	#include <wx/wx.h>
#else
	#include <wx/wxprec.h>
#endif

class mazeDlgApp : public wxApp
{
	public:
		bool OnInit();
		int OnExit();
};

#endif
