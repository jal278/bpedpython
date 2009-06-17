//---------------------------------------------------------------------------
//
// Name:        mazeApp.cpp
// Author:      Administrator
// Created:     2/5/2008 9:20:51 PM
// Description: 
//
//---------------------------------------------------------------------------

#include "mazeApp.h"
#include "mazeDlg.h"
#include "population.h"
#include <stdio.h>

IMPLEMENT_APP(mazeDlgApp)

bool mazeDlgApp::OnInit()
{
	mazeDlg* dialog = new mazeDlg(NULL);
	SetTopWindow(dialog);
	dialog->Show(true);		
	return true;
}
 
int mazeDlgApp::OnExit()
{
	return 0;
}
