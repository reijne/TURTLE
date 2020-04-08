//***********************************************************************
//
//	Name:			SharedMemory.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			20.05.1997
//
//
//
//
//
//***********************************************************************

#include "SharedMemory.h"

#include <stdlib.h>
#include <iostream>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <iomanip>
#include <unistd.h>
#include <time.h>

using namespace std;

const key_T InitialKey = 0x17930000 + geteuid() + (time(NULL) & 0x00FFFFFF);
key_T	SharedMemory::nextKey = InitialKey;

SharedMemory::SharedMemory(INT _size, Mode _mode)
{
	size = _size;
	tempSize = 0;
	mode = _mode;
	switch ( mode ) {
	case StandAlone:
		p = (char *) malloc(size);
		if (!p)
		{
                   cout << "Problem allocating memory" << endl;
		   cout << "Tried to allocate " << size << " bytes" << endl;
		   exit(1);
		}	

		break;
		
	case Master:
		key = nextKey++;
//		cout << "allocating shared memory..." << flush;
		id = shmget(key, size, IPC_CREAT | 0600);
//		cout << " key=" << hex << key << dec << ", id=" << id << "..." << flush;
		if ( id < 0 )
		{
			perror("SharedMemory::SharedMemory() master mode: shmget failed:");
			exit(1);
		}
		p = (char *) shmat(id, 0, 0);
//		cout << "p=" << hex << p << dec << endl;
		cout << endl;
		if ( p <= (void *) 0 )
		{
			perror("SharedMemory::SharedMemory() master mode: shmat failed:");
			exit(1);
		}
		freeShm();
		break;
		
	case Slave:
		break;
	}
}

SharedMemory::SharedMemory(key_T _key, INT _size)
{
	size = _size;
	tempSize = 0;
	key = _key;
	mode = Slave;
	id = shmget(key, size, 0);
	if ( id < 0 )
	{
		perror("SharedMemory::SharedMemory() slave mode: shmget failed:");
		exit(1);
	}
	p = (char *) shmat(id, 0, 0);
	if ( p <= (char *) 0 )
	{
		perror("SharedMemory::SharedMemory() slave mode: shmat failed:");
		exit(1);
	}
	freeShm();
}

void	SharedMemory::freeShm()
{
//	cout << "freeing shared memory: key=" << hex << key << dec
//		<< ", id=" << id << "..." << flush;
	if (shmctl(id, IPC_RMID, NULL)) 
	{
		perror("SharedMemory::freeShm() shmctl failed:");
		exit(1);
	}
//	cout << "o.k." << endl;
}


SharedMemory::~SharedMemory()
{
	switch ( mode ) {
	case StandAlone:
		free(p);
		break;
		
	case Master:
		if ( shmdt(p) == -1 )
		{
			perror("SharedMemory::~SharedMemory() master mode: shmdt failed:");
			exit(1);
		}
//		freeShm();
		break;
		
	case Slave:
		if ( shmdt(p) == -1 )
		{
			perror("SharedMemory::~SharedMemory() slave mode: shmdt failed:");
			exit(1);
		}
		break;
	}
}


void	SharedMemory::setSlave()
{
	mode = Slave;
}


void *	SharedMemory::allocate(INT n, INT _size)
{
char	*pp;

	if ( tempSize + n*_size > size )
	{
		cerr << "SharedMemory::allocate(): memory block exceeds reserved size" << endl;
		exit(1);
	}
	pp = p + tempSize;
	tempSize += n*_size;
	return ((void *) pp);
}




void	SharedMemory::deallocate(void *)
{
}



key_T	SharedMemory::getKey() const
{
	return key;
}

