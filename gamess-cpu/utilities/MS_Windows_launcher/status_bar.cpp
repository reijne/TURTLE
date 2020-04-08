#include <windows.h>
#include <windowsx.h>
#include <commctrl.h>
#include <string.h>
#include <stdio.h>
#include "status_bar.hpp"
#include "wgamess.hpp"
#include <stdlib.h>
#include "defines.h"


void status_bar::setstatus(status Stat,HWND hwnd)
{
    Status=Stat;
    showstatus(hwnd);
}

status status_bar::getstatus()
{
    return Status;
}

void status_bar::showstatus(HWND hwnd)
{
    char S[256];
   
    switch (Status)
    {
        case PROG_READY:
        case PROG_INIT:
        strcpy(S,"-- Ready --");
        break;

        case PROG_RUNNING:
        strcpy(S,"-- Running --");
        break;
    }    
    SendMessage(GetDlgItem(hwnd,ID_STATUSBAR),SB_SETTEXT,0,(LPARAM)S);
}

