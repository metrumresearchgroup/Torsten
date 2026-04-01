#include <gtest/gtest.h>
#include <stan/math/torsten/mpi/my_worker.hpp>
#include <iostream>
#include <vector>
#include <memory>

TEST(torsten_mpi_test, my_worker) {

  int np;
  int world_size = 3;
  int worker_id;

  // np > world size
  np = 10;
  // individual 0, 1, 2, 3 are worked on by rank 0
  for (int i = 0; i < 4; ++i) {
    worker_id = torsten::mpi::my_worker(i, np, world_size);
    EXPECT_EQ(worker_id, 0);
  }
  // individual 4, 5, 6, 7 are worked on by rank 1
  for (int i = 4; i < 7; ++i) {
    worker_id = torsten::mpi::my_worker(i, np, world_size);
    EXPECT_EQ(worker_id, 1);
  }

  // individual 7, 8, 9 are worked on by rank 2
  for (int i = 7; i < 10; ++i) {
    worker_id = torsten::mpi::my_worker(i, np, world_size);
    EXPECT_EQ(worker_id, 2);
  }

  // np < world size
  np = 2;
  worker_id = torsten::mpi::my_worker(0, np, world_size);
  EXPECT_EQ(worker_id, 0);
  worker_id = torsten::mpi::my_worker(1, np, world_size);
  EXPECT_EQ(worker_id, 1);

  // np = world size
  np = 3;
  worker_id = torsten::mpi::my_worker(0, np, world_size);
  EXPECT_EQ(worker_id, 0);
  worker_id = torsten::mpi::my_worker(1, np, world_size);
  EXPECT_EQ(worker_id, 1);
  worker_id = torsten::mpi::my_worker(2, np, world_size);
  EXPECT_EQ(worker_id, 2);
}
