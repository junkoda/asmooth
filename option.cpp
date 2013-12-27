#include <iostream>
#include <map>
#include <string>
#include <cstring>
#include "option.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;

Option::Option(int const argc, char const * const argv[])
{
  for(int j=1; j<argc-1; j+=2){
    if(argv[j][0] != '-')
      cout << "Error: " << argv[j] << " must be an option\n";
    else
      m_option.insert(pair<string,string>(argv[j], argv[j+1]));
  }
}

void Option::set_default(const char key[], const char value[])
{
  if(m_option.find(key) == m_option.end())
    m_option.insert(pair<string,string>(key, value));
}

void Option::set_default(const string key, const string value)
{
  if(m_option.find(key) == m_option.end())
    m_option.insert(pair<string,string>(key, value));
}
      

int Option::get_int(const char key[]) const
{
  return atoi(m_option.find(key)->second.c_str());
}

float Option::get_float(const char key[]) const
{
  return atof(m_option.find(key)->second.c_str());
}

void Option::print()
{
  for(map<string,string>::iterator p=m_option.begin(); p != m_option.end();
      ++p){
    cout << p->first << " " << p->second << endl;
  }
}
