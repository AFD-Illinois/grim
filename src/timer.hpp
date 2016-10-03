#include <sys/time.h>
#include <time.h>
#include <sstream>

// This file is mostly a copy from SpEC UtilsForTiming.*pp
// We define utilities allowing us to determine for how long the
// code has been running (useful for automated restarts)

class TimeStamp {
  friend double operator-(const TimeStamp &end,const TimeStamp &start);
public:
  TimeStamp(); // Initialized with current time
  // Initialized with number of seconds since 00:00 Jan 1 1970 GMT.
  double SecondsSinceEpoch() const;
private:
  long mSeconds,mMicroSeconds; // Since the epoch
};
double operator-(const TimeStamp &end,const TimeStamp &start);
const TimeStamp& TimeWhenExecutableStarted();

TimeStamp::TimeStamp() {
  timeval tv;
  gettimeofday(&tv,nullptr); // Second argument is a timezone ptr, unused/deprecated.
  mSeconds      = tv.tv_sec;
  mMicroSeconds = tv.tv_usec;
}

double operator-(const TimeStamp &end,const TimeStamp &start) {
  return static_cast<double>(end.mMicroSeconds-start.mMicroSeconds)/1000000.
    +    static_cast<double>(end.mSeconds-start.mSeconds);
}

double TimeStamp::SecondsSinceEpoch() const {
  return static_cast<double>(mMicroSeconds)/1000000. 
    +    static_cast<double>(mSeconds);
}

namespace {
  const TimeStamp TstartOfExecutable;
}

const TimeStamp& TimeWhenExecutableStarted() {
  return TstartOfExecutable;
}
