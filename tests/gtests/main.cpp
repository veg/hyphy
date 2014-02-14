#include "gtest/gtest.h"
#include <time.h>
#include "helperfunctions.h"

int main(int argc, char **argv)
{
  time_t timer;
  init_genrand (time(&timer));

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
