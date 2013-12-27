#ifndef _OPTION
#define _OPTION 1

#include <map>
#include <string>

class Option{
 public:
  Option(int const argc, char const * const argv[]);
  void set_default(const char key[], const char value[]);
  void set_default(const std::string key, const std::string value);
  int get_int(const char key[]) const;
  float get_float(const char key[]) const;
  void print(void);
  const std::string& operator [](const std::string& str) const {
    return m_option.find(str)->second;
  }
  const std::string& operator [](const char key[]) const {
    return m_option.find(key)->second;
  }
 private:
  std::map<std::string, std::string> m_option;
};

#endif
