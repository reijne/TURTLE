#include "../../../config.h"

#include <stdio.h>

#include "Fort31EndianConvert.h"


int	main(INT argc, char *argv[])
{
Fort31EndianConvert::Direction dir;
Fort31RecordFormat format;


	printf("program converts fort.31 file between little and big endian format\n");
	if ( argc!=5 )
	{	printf("usage: f31endian {l2b|b2l} {old|new|TRADPT} input-filename output-filename\n");
		return 1;
	}
	
	if ( !strcmp(argv[1], "b2l") )
		dir = Fort31EndianConvert::BigToLittle;
	else
	if ( !strcmp(argv[1], "l2b") )
		dir = Fort31EndianConvert::LittleToBig;
	else
	{
		printf("unknown direction \"%s\"\n", argv[1]);	
		exit(1);
	}
	
	if ( !strcmp(argv[2], "old") )
		format = Fort31RecordFormatOld;
	else
	if ( !strcmp(argv[2], "new") )
		format = Fort31RecordFormatNew;
	else
	if ( !strcmp(argv[2], "TRADPT") )
		format = Fort31RecordFormatTRADPT;
	else
	{
		printf("unknown format \"%s\"\n", argv[2]);	
		exit(1);
	}
	
	
	
	printf("converting...\n");
Fort31EndianConvert	fcnv(argv[3], argv[4], dir, format);	
	printf("ready.\n");
}
