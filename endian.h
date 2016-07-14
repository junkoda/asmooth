// endian_int() converts an int read from file to correct endian
// same for short and float


#ifdef SWAPENDIAN

union Int 
{
  char b[4];
  int num;
};

union Short
{
  char b[2];
  short num;
};

union Float
{
  char b[4];
  float num;
};

static inline int endian_int(const int x)
{
  Int temp;
  char const * const p= (char const*) &x;
  temp.b[0]= p[3];
  temp.b[1]= p[2];
  temp.b[2]= p[1];
  temp.b[3]= p[0];

  return temp.num;
}

static inline short endian_short(const short x)
{
  Short temp;
  char const * const p= (char const*) &x;
  temp.b[0]= p[1];
  temp.b[1]= p[0];

  return temp.num;
}

static inline float endian_float(const float x)
{
  Float temp;
  char const * const p= (char const *) &x;
  temp.b[0]= p[3];
  temp.b[1]= p[2];
  temp.b[2]= p[1];
  temp.b[3]= p[0];

  return temp.num;
}

#else
static inline int endian_int(const int x) { return x; }
static inline short endian_short(const short x) { return x; }
static inline float endian_float(const float x) { return x; }
#endif

