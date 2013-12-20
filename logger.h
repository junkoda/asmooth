#ifndef LOGGER_H
#define LOGGER_H 1

#include <ostream>
#include <sstream>
#include "performance_items.h"

class Logger : public std::ostringstream {
 public:
  Logger();
  //~PerformanceLog();
  void begin_timer(performance_item item);
  void end_timer();
  void print_log(std::ostream& out);
  //std::ostream& operator<<(std::ostream& left);
 private:
  double times[performance_item(last)];
  double last_time;
  performance_item prev_item;
  //ostringstream log_buf;
  
  static double now();
};

#endif
