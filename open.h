#ifndef OPEN_H
#define OPEN_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/stat.h>

//
// Open file like fopen, but create directory if necesary.
//
FILE* open(const char filename[], const char mode[]);
int make_directory(const char path[]);

#endif
