#include <direct.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <process.h>
#include <errno.h>

static char dfltin[]=".in";
static char dfltout[]=".out";

void main(int argc, char *argv[ ], char *envp[ ] ) {
	/* Program to replace the UNIX shell script "rungamess" on Windows machines.
	   Originally a batch script was attempted but the limitations in the language
	   and the portability issues (every version of the command interpreter is slightly
	   different) rendered this approach completely hopeless.

	   The program does needs to construct an input file name and an output file name. 
	   We need to check that a scratch directory has been defined and exists.
	   We need a jobname to build the final scratch directory and create that.
	   We need to run GAMESS-UK in the scratch directory using the input and output files.
	   We need to tidy up and remove the scratch directory.

	   That's all folks!
   */
	FILE *tmp;
    char cwd[_MAX_PATH];
	char inputfile[_MAX_PATH];
	char outputfile[_MAX_PATH];
	char jobname[_MAX_PATH];
	char scratchdir[_MAX_PATH];
	char gamessbin[_MAX_PATH];
	char gamessexe[_MAX_PATH];
	size_t jb_len;
	size_t fl_len;
	char *cptr;
	char command[3*_MAX_PATH];
	intptr_t exit_state;
	int exit_error;

	if( _getcwd( cwd, _MAX_PATH ) == NULL ) {
      perror( "_getcwd error" );
	  exit(5);
	};
    
	if (!(argc==2 || argc==3)) {
		printf("Usage: rungamess <inputfile> [<outputfile>]\n\n");
		printf("Runs GAMESS-UK using <inputfile> and writing <outputfile>.\n\n");
		printf("It is assumed that <inputfile> ends in .in\n");
		printf("if this is not the case <inputfile>.in will be used for input.\n\n");
		printf("It is assumed that <outputfile> ends in .out\n");
		printf("if this is not the case <outputfile>.out will be used for output.\n");
		printf("If no output file was specified the output file name will be \n");
		printf("constructed from the input file name\n");
		exit(10);
	}

	// Sort out input file name

	if (argv[1][1]==':') {
		strcpy(inputfile,argv[1]);
	} else {
		strcpy(inputfile,cwd);
		strcat(inputfile,"\\");
		strcat(inputfile,argv[1]);
	}
	fl_len = strlen(inputfile);
	if (strcmp(&inputfile[fl_len-sizeof(dfltin)+1],dfltin)!=0) {
		strcat(inputfile,dfltin);
		fl_len = strlen(inputfile);
	};
	if ((tmp = fopen(inputfile,"r"))==NULL){
		printf("Input file: %s\n",inputfile);
		printf("does not exist or cannot be read.\n");
		exit(20);
	};
	fclose(tmp);

	// Sort out output file name
	
	if (argc==3) {
	    if (argv[2][1]==':') {
		    strcpy(outputfile,argv[2]);
	    } else {
		    strcpy(outputfile,cwd);
			strcat(outputfile,"\\");
		    strcat(outputfile,argv[2]);
	    } 
	    fl_len = strlen(outputfile);
	    if (strcmp(&outputfile[fl_len-sizeof(dfltout)+1],dfltout)!=0) {
		    strcat(outputfile,dfltout);
	    }
	} else {
		strcpy(outputfile,inputfile);
		strcpy(&outputfile[fl_len-sizeof(dfltin)+1],dfltout);
	};
	if ((tmp = fopen(outputfile,"r+"))==NULL){
		if ((tmp = fopen(outputfile,"w"))==NULL){
		    printf("Output file: %s\n",outputfile);
		    printf("cannot be created or is not writable.\n");
		    exit(30);
		};
	};
	fclose(tmp);

	// Sort out the jobname

	cptr = strrchr(inputfile,'\\');
	strcpy(jobname,++cptr);
	jb_len = strlen(jobname);
	jobname[jb_len-sizeof(dfltin)+1]='\0';

	// Sort out the name of the bin directory

    if ((cptr = getenv("GAMESS_BIN"))==0) {
		printf("Environment variable GAMESS_BIN not defined\n");
		printf("This variable is needed to construct the executable name!\n");
		exit(40);
	}
	strcpy(gamessbin,cptr);
	if (_chdir(gamessbin)!=0) {
		printf("The directory %s\n",gamessbin);
		printf("defined in environment variable GAMESS_BIN does not exist!\n");
		exit(50);
	};

	// Sort out the name of the scratch directory

	if ((cptr = getenv("GAMESS_SCR"))==0) {
		printf("Environment variable GAMESS_SCR not defined\n");
		printf("This variable is needed to construct the scratch directory!\n");
		exit(60);
	}
	strcpy(scratchdir,cptr);
	if (_chdir(scratchdir)!=0) {
		printf("The directory %s\n",scratchdir);
		printf("defined in environment variable GAMESS_SCR does not exist!\n");
		exit(70);
	};
	strcat(scratchdir,jobname);
	_mkdir(scratchdir);
	if (_chdir(scratchdir)!=0) {
		printf("Could not create directory %s\n",scratchdir);
		printf("Perhaps there is a file with the same name or\n");
		printf("you do not have write access permissions.\n");
		exit(80);
	};

	// Build and run the command

	sprintf(gamessexe,"GAMESS_EXE=\"%s\\gamess.exe\"",gamessbin);
	_putenv(gamessexe);
	sprintf(command,"%%GAMESS_EXE%% < \"%s\" > \"%s\"",inputfile,outputfile);
	//printf("%s\n",gamessexe);
	//printf("Running: %s\n",command);
	exit_state = system(command);
	exit_error = errno;
	if (exit_state == 0) { }
	else if ((exit_state == -1) && (exit_error == E2BIG)) {
		printf("The argument list is too big\nThe command was: %s",command);
	}
	else if ((exit_state == -1) && (exit_error == ENOENT)) {
		printf("The command interpreter cmd.exe was not found.");
	}
	else if ((exit_state == -1) && (exit_error == ENOEXEC)) {
		printf("The command interpreter file cmd.exe has an invalid format and is not executable.");
	}
	else if ((exit_state == -1) && (exit_error == ENOMEM)) {
		printf("There is not enough memory available to run the command\nor the memory is corrupted.");
	}
	else if (exit_state == 1) {
		printf("The command: %s was not found\nCould not start the job.",gamessexe);
	}
	else {
		printf("GAMESS-UK job failed\n");
		printf("exit code %d\n",exit_state);
		printf("error number %d\n",errno);
		exit(exit_state);
	};

	// Tidy up the scratch directory

	_chdir(cptr);
    sprintf(command,"DEL /Q \"%s\"",scratchdir);
	exit_state = system(command);
    if (exit_state != 0) {
		printf("Failed to delete scratch files\n");
		printf("exit code %d\n",exit_state);
		printf("error number %d\n",errno);
		if (exit_state > 0) {
			exit(exit_state);
		} else {
			exit(errno);
		};
	};
	if (_rmdir(scratchdir)!=0) {
		printf("Failed to delete scratch directory %s\n",scratchdir);
		perror("_rmdir error");
		exit(90);
	};

}
