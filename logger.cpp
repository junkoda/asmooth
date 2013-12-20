#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include "logger.h"

using namespace std;

Logger::Logger()
{
  const int n= performance_item(last);
  for(int j=0; j<n; ++j)
    times[j]= 0.0;

  prev_item= etc;
  last_time= now();
}

void Logger::begin_timer(performance_item item)
{
  const double time_now= now();

  times[prev_item] += (time_now - last_time);
  last_time= time_now;
  prev_item= item;
}

void Logger::end_timer()
{
  begin_timer(performance_item(etc));
};

void Logger::print_log(ostream& out)
{
  //out << ostringstream::str() << endl;
  out << ostringstream::str();
  ostringstream::str("");

  //out.setf(ios_base::scientific, ios_base::floatfield);
  //out.setf(ios_base::left, ios_base::adjustfield);

  const int n= performance_item(last);
  double total= 0;
  for(int j=0; j<n; ++j)
    total += times[j];

  char buf[128];
  for(int j=0; j<n; ++j) {
    /*
    out << "time "
        << ios_base::setw(20) << left << performance_name[j] << " "
        << setprecision(6) << scientific 
        << times[j] << "  "
        << setw(6) << fixed << setprecision(2) << right 
	<< times[j]/total*100 << "%" << endl;
    */
    sprintf(buf, "time %20s %8.2e %4.1f%%\n",
	    performance_name[j], times[j], times[j]/total*100);
    times[j]= 0.0;

    out << buf;
  }
  sprintf(buf, "time %20s %10.4e\n",
	  "total", total);
  out << buf;

}

double Logger::now()
{
  struct timeval tp;
  struct timezone tzp;
  const int i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
