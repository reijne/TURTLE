void GAMESS_OnDestroy(HWND hwnd);
BOOL GAMESS_OnCreate(HWND hwnd,CREATESTRUCT FAR * lpCreateStruct);
BOOL Register(HINSTANCE hInst);
HWND Create(HINSTANCE hInst,int nCmdShow);
LRESULT CALLBACK WndProc(HWND hwnd,UINT Message,WPARAM wParam,LPARAM lParam);
void GAMESS_OnCommand(HWND hwnd,int id,HWND hwndCtl,UINT codeNotify);

#define HANDLE_MM_MCINOTIFY(hwnd,wParam,lParam,fn)\
((fn)((hwnd),(UINT)(wParam),(int)LOWORD(lParam)),0L)

void OpenFileName(HWND hwnd,char S[256],int *ret);
void SaveFileName(HWND hwnd,char S[256],int *ret);
extern "C" gamess(char input_file[],char output_file[]);
void rungamess();
DWORD gogamess();

void getdate(char date[], char time[], char year[]);
