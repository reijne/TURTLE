/* $Header: /c/qcg/cvs/psh/GAMESS-UK/g/tcgmsg/ipcv4.0/sockets.h,v 1.1.1.6 2007-10-30 10:14:16 jmht Exp $ */

extern void ShutdownAll();
extern int ReadFromSocket();
extern int WriteToSocket();
extern void CreateSocketAndBind();
extern int ListenAndAccept();
extern int CreateSocketAndConnect();
extern long PollSocket();
extern long WaitForSockets(int nsock, int *socks, int *list);
