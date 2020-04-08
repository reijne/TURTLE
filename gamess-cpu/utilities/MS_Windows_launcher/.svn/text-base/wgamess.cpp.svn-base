#include <windows.h>
#include <windowsx.h>
#include <winbase.h>
// #include <env.h>
#include <stdlib.h>
#include <commctrl.h>
#include <iostream.h>
#include <stdio.h>
#include <time.h>
#include <process.h>
#include <errno.h>
#include "wgamess.hpp"
#include "status_bar.hpp"
#include "defines.h"
#include "resource.h"

HWND hwnd;
static HINSTANCE hInstance;
static char szAppName[]="Gamess";
status_bar program_state;
HANDLE Gamess_thread;

static char defaultin[]=".in";
static char defaultout[]=".out";

#define LEN_FILENAME 512
char input_file[LEN_FILENAME],output_file[LEN_FILENAME];
char jobname[_MAX_PATH];
char gamessrun[_MAX_PATH];
#define NO_FILENAMES 162
char *filenames[]={
"ed0","ed1","ed2" ,"ed3" ,"ed4" ,"ed5" ,"ed6" ,
"ed7" ,"ed8" ,"ed9" ,"ed10","ed11","ed12","ed13",
"ed14","ed15","ed16","ed17","ed18","ed19",
"mt0" ,"mt1" ,"mt2" ,"mt3" ,"mt4" ,"mt5" ,"mt6" ,
"mt7" ,"mt8" ,"mt9" ,"mt10","mt11","mt12","mt13",
"mt14","mt15","mt16","mt17","mt18","mt19",       
"ft01","ft02","ft03","ft04","ft05","ft06","ft07","ft08",
"ft09","ft10","ft11",
"ft12","ft13","ft14","ft15","ft16","ft17","ft18","ft19",
"ft20","ft21","ft22","ft23","ft24","ft25","ft26","ft27",
"ft28","ft29","ft30","ft31",
"ft32","ft33","ft34","ft35","ft36","ft37","ft38","ft39",
"ft40","ft41","ft42","ft43","ft44","ft45","ft46","ft47",
"ft48","ft49","ft50","ft51",
"ft52","ft53","ft54","ft55","ft56","ft57","ft58","ft59",
"ft60",
"ftn001","ftn002","ftn003","ftn004","ftn005","ftn006","ftn007",
"ftn008","ftn009","ftn010","ftn011",
"ftn012","ftn013","ftn014","ftn015","ftn016","ftn017","ftn018",
"ftn019","ftn020","ftn021","ftn022","ftn023","ftn024","ftn025",
"ftn026","ftn027","ftn028","ftn029","ftn030","ftn031",
"ftn032","ftn033","ftn034","ftn035","ftn036","ftn037","ftn038",
"ftn039",
"ftn040","ftn041","ftn042","ftn043","ftn044","ftn045",
"ftn046","ftn047","ftn048","ftn049","ftn050","ftn051",
"ftn052","ftn053","ftn054","ftn055","ftn056","ftn057",
"ftn058","ftn059","ftn060","table","table-ci"
};

#define NO_ENV 165
char *genv[NO_ENV];
#define NO_ARG 4
char *gargv[NO_ARG];

extern int errno;

int __stdcall WinMain(HINSTANCE hInst, HINSTANCE hPrev, LPSTR lpszCmdParam,int nCmdShow)
{
    MSG msg;

    if (!hPrev)
       if (!Register(hInst)) return FALSE;

    hwnd=Create(hInst,nCmdShow);
    if (!hwnd) return FALSE;

    while (GetMessage(&msg,NULL,0,0))
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    return msg.wParam;
}

BOOL Register(HINSTANCE hInst)
{
    WNDCLASS wnd;    

    wnd.style=CS_HREDRAW|CS_VREDRAW;
    wnd.lpfnWndProc=WndProc;
    wnd.cbClsExtra=0;
    wnd.cbWndExtra=DLGWINDOWEXTRA;
    wnd.hInstance=hInst;
    wnd.hIcon=LoadIcon(NULL,IDI_APPLICATION);
    wnd.hCursor=LoadCursor(NULL,IDC_ARROW);
    wnd.hbrBackground=GetStockBrush(LTGRAY_BRUSH);
    wnd.lpszMenuName=NULL;
    wnd.lpszClassName=szAppName;

    return (RegisterClass(&wnd)!=0);
}

HWND Create(HINSTANCE hInst,int nCmdShow)
{
    HWND hwnd;
    
    hInstance=hInst;
    
    hwnd=CreateDialog(hInst,szAppName,0,NULL);

    if (hwnd==NULL)
       return FALSE;

    ShowWindow(hwnd,nCmdShow);        
    program_state.setstatus(PROG_INIT,hwnd); 
	if (program_state.getstatus()==PROG_INIT)
    {
        program_state.setstatus(PROG_READY,hwnd);
        for (int i=0;i<NO_FILENAMES;i++)ComboBox_AddString(GetDlgItem(hwnd,ID_COMBO_ENV),filenames[i]);          
    }
    return hwnd;
}

LRESULT CALLBACK WndProc(HWND hwnd,UINT Message,WPARAM wParam,LPARAM lParam)
{
    switch (Message)
    {
        case WM_CTLCOLORDLG:
        return (BOOL)GetStockBrush(LTGRAY_BRUSH);
        break;

        case WM_CTLCOLORBTN:
        case WM_CTLCOLORSTATIC:
        SetTextColor((HDC)wParam,RGB(0,0,0));
        SetBkMode((HDC)wParam,TRANSPARENT);
        return(BOOL)GetStockBrush(LTGRAY_BRUSH);     
        break;
        
        HANDLE_MSG(hwnd,WM_CREATE,GAMESS_OnCreate);
        HANDLE_MSG(hwnd,WM_DESTROY,GAMESS_OnDestroy);
        HANDLE_MSG(hwnd,WM_COMMAND,GAMESS_OnCommand);
        default:
            return DefWindowProc(hwnd,Message,wParam,lParam);
    }
}

BOOL GAMESS_OnCreate(HWND hwnd,CREATESTRUCT FAR * lpCreateStruct)
{    
    InitCommonControls();  
    return TRUE;
}

void GAMESS_OnDestroy(HWND hwnd)
{
    PostQuitMessage(0);
}

void GAMESS_OnCommand(HWND hwnd,int id,HWND hwndCtl,UINT codeNotify)
{
    int ret;
    char S[LEN_FILENAME];
    char internal[LEN_FILENAME],external[LEN_FILENAME];
    char date[6],time[9],year[5];
    char date_time[30];
    int *index_list;
    int no_items;   

    switch (id)
    {        
        case ID_SEL_INPUT:
        if (program_state.getstatus()==PROG_READY||program_state.getstatus()==PROG_INIT)
        {
            OpenFileName(hwnd,input_file,&ret);
            if (ret==0)
            {
				size_t fl_len;
				size_t in_len = strlen(defaultin);
				size_t out_len = strlen(defaultout);
                strcpy(output_file,input_file);
				fl_len = strlen(output_file);
                if (strcmp(&output_file[fl_len-in_len],defaultin)==0)
				{
					strcpy(&output_file[fl_len-in_len],defaultout);
				} else
				{
                    strcat(output_file,defaultout);
				};
				SetWindowText(GetDlgItem(hwnd,IDC_EDIT1),input_file);
                SetWindowText(GetDlgItem(hwnd,IDC_EDIT2),output_file);
            }
        } else
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };            
        break;
        
        case ID_SEL_OUTPUT:
        if (program_state.getstatus()==PROG_READY||program_state.getstatus()==PROG_INIT)
        {
            SaveFileName(hwnd,output_file,&ret);
            if (ret==0)
            {
                SetWindowText(GetDlgItem(hwnd,IDC_EDIT2),output_file);
            }
        }
        else
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };              
        break;
        
		case IDC_EDIT1:
        if (program_state.getstatus()==PROG_READY||program_state.getstatus()==PROG_INIT)
        {
			GetWindowText(GetDlgItem(hwnd,IDC_EDIT1),input_file,LEN_FILENAME);
			size_t fl_len;
			size_t in_len = strlen(defaultin);
			size_t out_len = strlen(defaultout);
            strcpy(output_file,input_file);
			fl_len = strlen(output_file);
            if (strcmp(&output_file[fl_len-in_len],defaultin)==0)
			{
				strcpy(&output_file[fl_len-in_len],defaultout);
			} else
			{
                strcat(output_file,defaultout);
			};
            SetWindowText(GetDlgItem(hwnd,IDC_EDIT2),output_file);
        } else
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };            
        break;
        
        case IDC_EDIT2:
        if (program_state.getstatus()==PROG_READY||program_state.getstatus()==PROG_INIT)
        {
            GetWindowText(GetDlgItem(hwnd,IDC_EDIT2),output_file,LEN_FILENAME);
        }
        else
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };              
        break;

        case ID_START:
        if (program_state.getstatus()==PROG_READY||program_state.getstatus()==PROG_INIT)
        {
            if (strlen(input_file)!=0&&strlen(output_file)!=0)
            {
                // program_state.setstatus(PROG_RUNNING,hwnd);
                getdate(date,time,year);
                sprintf(date_time,"%s-%s at %s",date,year,time);
                SetWindowText(GetDlgItem(hwnd,ID_TEXT_DATE),date_time);
				// DEBUG
				//MessageBox(hwnd,input_file,"Information",MB_OK);
				//MessageBox(hwnd,output_file,"Information",MB_OK);
				// DEBUG
                rungamess();
            }else
            {
                MessageBox(hwnd,"No input or output file supplied","Error",MB_OK|MB_ICONSTOP);
            };            
        }
        else
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };          
        break;
       
        case ID_SAVE_EXT:
        if (program_state.getstatus()==PROG_READY||program_state.getstatus()==PROG_INIT)
        {
            SaveFileName(hwnd,S,&ret);
            if (ret==0)
            {
				SetWindowText(GetDlgItem(hwnd,IDC_EDIT3),S);
            }         
        } else
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };  
        break;
        
		case IDC_EDIT3:
        if (program_state.getstatus()==PROG_READY||program_state.getstatus()==PROG_INIT)
        {
			GetWindowText(GetDlgItem(hwnd,IDC_EDIT3),S,LEN_FILENAME);
                    
        } else
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };  
        break;

        case ID_COMBO_ENV:
        if (program_state.getstatus()==PROG_INIT || program_state.getstatus()==PROG_READY)
        {
            if (codeNotify==CBN_SELCHANGE)
            {
                ComboBox_GetLBText(GetDlgItem(hwnd,ID_COMBO_ENV),ComboBox_GetCurSel(GetDlgItem(hwnd,ID_COMBO_ENV)),internal);
                SetWindowText(GetDlgItem(hwnd,IDC_EDIT3),internal);
            }
        } else
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };  
        break;
        
        case ID_LIST_ENV:
        if (program_state.getstatus()==PROG_READY||program_state.getstatus()==PROG_INIT)
        {            
        } else
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };  
        break;
        
        case ID_CLEAR_ENV:
        if (program_state.getstatus()==PROG_READY||program_state.getstatus()==PROG_INIT)
        {
            no_items=ListBox_GetSelCount(GetDlgItem(hwnd,ID_LIST_ENV));
            for (;no_items>0;){                
                index_list=new int [no_items];
                ListBox_GetSelItems(GetDlgItem(hwnd,ID_LIST_ENV), no_items, index_list);          
                ListBox_GetText(GetDlgItem(hwnd,ID_LIST_ENV), *index_list, S);
                sscanf(S,"%s",internal);
                // setenv(internal,NULL,1);
                ListBox_DeleteString(GetDlgItem(hwnd,ID_LIST_ENV), *index_list);
                delete [] index_list;
                no_items=ListBox_GetSelCount(GetDlgItem(hwnd,ID_LIST_ENV));
            };            
        } else  
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };  
        break;
        
        case ID_SET_ENV:
        if (program_state.getstatus()==PROG_READY)
        {            
            if (ComboBox_GetCurSel(GetDlgItem(hwnd,ID_COMBO_ENV))>-1) {
                ComboBox_GetLBText(GetDlgItem(hwnd,ID_COMBO_ENV),ComboBox_GetCurSel(GetDlgItem(hwnd,ID_COMBO_ENV)),internal);
                GetWindowText(GetDlgItem(hwnd,IDC_EDIT3),external,sizeof(external));
                // setenv(internal,external,1);
                sprintf(S,"%s=%s",internal,external);
                ListBox_AddString(GetDlgItem(hwnd,ID_LIST_ENV),S);
            } else
            {
                MessageBox(hwnd,"Please first select internal filename","Error",MB_OK|MB_ICONSTOP);
            };
                
        } else if (program_state.getstatus()==PROG_INIT)
        {
            MessageBox(hwnd,"Please first select internal filename","Error",MB_OK|MB_ICONSTOP);
        } else
        {
            MessageBox(hwnd,"You cannot do this while a job is running","Error",MB_OK|MB_ICONSTOP);
        };    
        break;      
    }
}

void OpenFileName(HWND hwnd,char S[256],int *ret)
{
    OPENFILENAME ofn;
    char szDirName[256];
    char szFileTitle[256];    
    char szFile[256];       
   
    szFile[0]=0;
    memset(&ofn,0,sizeof(OPENFILENAME));
    ofn.lStructSize=sizeof(OPENFILENAME);
    ofn.hwndOwner=hwnd;    
    ofn.lpstrFilter="All files (*.*) \0*.*\0";       
    ofn.nFilterIndex=1;
    ofn.lpstrFile=szFile;
    ofn.nMaxFile=256;
    ofn.lpstrFileTitle=szFileTitle;
    ofn.nMaxFileTitle=sizeof(szFileTitle);
    ofn.lpstrInitialDir=szDirName;
    ofn.Flags=OFN_FILEMUSTEXIST;
    *ret=0;
    if (GetOpenFileName(&ofn)) strcpy(S,ofn.lpstrFile);
    else *ret=1;            
}

void SaveFileName(HWND hwnd,char S[256],int *ret)
{
    OPENFILENAME ofn;
    char szDirName[256];
    char szFile[256],szFileTitle[256];        

    memset(&ofn,0,sizeof(OPENFILENAME));
    szFile[0]='\0';
    *ret=0;
    ofn.lStructSize=sizeof(OPENFILENAME);
    ofn.hwndOwner=hwnd;
    ofn.lpstrFilter="All files (*.*)\0*.*\0";
    ofn.lpstrFile=szFile;
    ofn.nMaxFile=sizeof(szFile);
    ofn.lpstrFileTitle=szFileTitle;
    ofn.nMaxFileTitle=sizeof(szFileTitle);
    ofn.lpstrInitialDir=szDirName;

    ofn.Flags=OFN_OVERWRITEPROMPT;

    if (GetSaveFileName(&ofn)) strcpy(S,ofn.lpstrFile);
    else *ret=1;
}

void rungamess()
{
    DWORD TID;
    Gamess_thread=CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)gogamess,NULL,0,&TID);
    SetThreadPriority(Gamess_thread,THREAD_PRIORITY_HIGHEST);   
}

DWORD gogamess()
{
	intptr_t exit_state;
	int exit_error;
	char *cptr;
	char command[3*LEN_FILENAME];

    // gamess(input_file,output_file); 
	program_state.setstatus(PROG_RUNNING,hwnd);

	// setup the command line arguments
	cptr = strrchr(input_file,'\\');
	strcpy(jobname,++cptr);


	// setup the environment variables
	//_putenv("GAMESS_GUI=active");
	int no_items = ListBox_GetCount(GetDlgItem(hwnd,ID_LIST_ENV));
	//genv[no_items+1]=NULL;
	//genv[no_items]=new char [LEN_FILENAME];
	//strcpy(genv[no_items],"GAMESS_GUI=active");
	while (no_items>0) {
		char envvar[LEN_FILENAME];
		//genv[no_items-1] = new char [LEN_FILENAME];
		ListBox_GetText(GetDlgItem(hwnd,ID_LIST_ENV),no_items-1,envvar);
		_putenv(envvar);
		no_items--;
	}

	//exit_state = _spawnvpe(_P_WAIT,"rungamess.bat",gargv,genv);

	// build and run the command
	cptr = getenv("GAMESS_BIN");
	if (cptr == NULL) {
		exit_state = 40;
	} 
	else {
	    sprintf(gamessrun,"GAMESS_RUN=\"%s\\rungamess.exe\"",cptr);
	    _putenv(gamessrun);
		sprintf(command,"%%GAMESS_RUN%% \"%s\" \"%s\"",input_file,output_file);
	    // DEBUG
	    //MessageBox(hwnd,gamessrun,"Information",MB_OK);
	    //MessageBox(hwnd,command,"Information",MB_OK);
	    // DEBUG
        exit_state = system(command);
		exit_error = errno;

    // The spawn call never worked as it kept complaining about _P_WAIT having an illegal
	// value...
	//	sprintf(command,"\"%s\\rungamess.exe\"",cptr);
	//    // DEBUG
	//    MessageBox(hwnd,command,"Information",MB_OK);
	//    // DEBUG _P_WAIT
    //    exit_state = _spawnl(_P_WAIT,command,command,input_file,output_file,NULL);
	//	exit_error = errno;
	}

	// clear the environment variables
	no_items = ListBox_GetCount(GetDlgItem(hwnd,ID_LIST_ENV));
	while (no_items>0) {
		char envvar[LEN_FILENAME];
		ListBox_GetText(GetDlgItem(hwnd,ID_LIST_ENV),no_items-1,envvar);
		cptr = strrchr(envvar,'=');
		*(++cptr)='\0';
		_putenv(envvar);
		no_items--;
	}

	program_state.setstatus(PROG_READY,hwnd);
	if (exit_state==0) {
        MessageBox(hwnd,"GAMESS job finished","Information",MB_OK);
	}
	else if ((exit_state==-1)&&(exit_error==E2BIG)) {
		char msg[512];
		sprintf(msg,"The argument list is too big\nThe command was: %s",gamessrun);
		MessageBox(hwnd,msg,"Error",MB_OK|MB_ICONSTOP);
	}
	else if ((exit_state==-1)&&(exit_error==ENOENT)) {
		MessageBox(hwnd,"The command interpreter cmd.exe was not found.","Error",MB_OK|MB_ICONSTOP);
	}
	else if ((exit_state==-1)&&(exit_error==ENOEXEC)) {
		MessageBox(hwnd,"The command interpreter file cmd.exe has an invalid format and is not executable.","Error",MB_OK|MB_ICONSTOP);
	}
	else if ((exit_state==-1)&&(exit_error==ENOMEM)) {
		MessageBox(hwnd,"There is not enough memory available to run the command\nor the memory is corrupted.","Error",MB_OK|MB_ICONSTOP);
	}
	else if (exit_state==1) {
        char msg[512];
		sprintf(msg,"The command: %s was not found\nCould not start the job.",gamessrun);
		MessageBox(hwnd,msg,"Error",MB_OK|MB_ICONSTOP);
	} 
	else if (exit_state==5) {
		MessageBox(hwnd,"rungamess.exe could not obtain the current working directory",
			       "Error",MB_OK|MB_ICONSTOP);
	}
	else if (exit_state==10) {
		char msg[512];
		sprintf(msg,"The command: %s\ndoes not conform to the rungamess.exe usage.",command);
		MessageBox(hwnd,msg,"Error",MB_OK|MB_ICONSTOP);
	}
	else if (exit_state==20) {
		char msg[512];
		sprintf(msg,"The input_file: %s\ncould not be read.",input_file);
		MessageBox(hwnd,msg,"Error",MB_OK|MB_ICONSTOP);
	}
	else if (exit_state==30) {
		char msg[512];
		sprintf(msg,"The output_file: %s\ncould not be written.",output_file);
		MessageBox(hwnd,msg,"Error",MB_OK|MB_ICONSTOP);
	}
	else if (exit_state==40) {
		MessageBox(hwnd,"The environment variable GAMESS_BIN appears not to have been set.\nPlease check your environment settings.","Error",MB_OK|MB_ICONSTOP);
	}
	else if (exit_state==50) {
		MessageBox(hwnd,"The directory defined by environment variable GAMESS_BIN cannot be accessed.","Error",MB_OK|MB_ICONSTOP);
	}
	else if (exit_state==60) {
		MessageBox(hwnd,"The environment variable GAMESS_SCR (scratch directory) appears not to have been set.\nPlease check your environment settings.","Error",MB_OK|MB_ICONSTOP);
	}
	else if (exit_state==70) {
		MessageBox(hwnd,"The directory defined by environment variable GAMESS_SCR (scratch directory) cannot be accessed.","Error",MB_OK|MB_ICONSTOP);
	}
	else if (exit_state==80) {
		MessageBox(hwnd,"Could not create a scratch sub-directory in the directory defined by GAMESS_SCR.","Error",MB_OK|MB_ICONSTOP);
	}
	else if (exit_state==90) {
		MessageBox(hwnd,"Could not delete the scratch sub-directory in the directory defined by GAMESS_SCR after the job finished.","Error",MB_OK|MB_ICONSTOP);
	}
	else {
		char msg[512];
		sprintf(msg,"Error detected!\nGamess job finished\nExit code: %d\nError code: %d",exit_state,errno);
		MessageBox(hwnd,msg,"Error",MB_OK|MB_ICONSTOP);
	}
    // SendMessage(hwnd,WM_DESTROY,0,0);    
    ExitThread(1);
    return 0;
}

extern "C" void winexit()
{        
    MessageBox(hwnd,"Error detected!\nGamess job finished","Error",MB_OK|MB_ICONSTOP);
    // SendMessage(hwnd,WM_DESTROY,0,0);    
    TerminateThread(Gamess_thread,1);
}

void getdate(char zdate[], char ztime[], char zyear[])
{
	time_t my_time;
	struct tm tp;

	time(&my_time);
	tp = *localtime(&my_time);

	sprintf(ztime,"%02i:%02i:%02i",tp.tm_hour,tp.tm_min,tp.tm_sec);
	sprintf(zdate,"%02i-%02i",tp.tm_mday,tp.tm_mon+1);
	sprintf(zyear,"%i",tp.tm_year+1900);
}
