#include "adamic_adar_test.h"
#include "ctime"
#define restrict

class AdamicAdarTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    stinger_config = (struct stinger_config_t *)xcalloc(1,sizeof(struct stinger_config_t));
    stinger_config->nv = 1<<13;
    stinger_config->nebs = 1<<16;
    stinger_config->netypes = 2;
    stinger_config->nvtypes = 2;
    stinger_config->memory_size = 1<<30;
    S = stinger_new_full(stinger_config);
    xfree(stinger_config);

    stinger_insert_edge_pair(S,0,0,1,1,1);
    stinger_insert_edge_pair(S,0,0,2,1,1);
    stinger_insert_edge_pair(S,0,0,3,1,1);
    stinger_insert_edge_pair(S,0,1,4,1,1);
    stinger_insert_edge_pair(S,0,2,4,1,1);
    stinger_insert_edge_pair(S,0,1,5,1,1);
    stinger_insert_edge_pair(S,0,1,6,1,1);
    stinger_insert_edge_pair(S,0,1,7,1,1);
    stinger_insert_edge_pair(S,0,3,8,1,1);
  }

  virtual void TearDown() {
    stinger_free_all(S);
  }

  struct stinger_config_t * stinger_config;
  struct stinger * S;
};

TEST_F(AdamicAdarTest, AllEdgeTypes) {
  int64_t * candidates = NULL;
  double * scores = NULL;
  int64_t * candidates2 = NULL;
  double * scores2 = NULL;
  int64_t len = adamic_adar(S, 0, -1, &candidates, &scores);
  printf("\n\n\n\nwhy \n\n\n\n");
  int64_t len2 = adamic_adar(S, 4, -1, &candidates2, &scores2);
  for (int64_t i = 0; i < len; i++) {
    switch (candidates[i]) {
      case 0:
        EXPECT_DOUBLE_EQ(1.0/log(2)+1/log(5),scores[i]);
        break;
      case 1:
        EXPECT_DOUBLE_EQ(1.0/log(2)+1/log(5),scores[i]);
        break;
      case 2:
        EXPECT_DOUBLE_EQ(1.0/log(2)+1/log(5),scores[i]);
        break;
      case 3:
        EXPECT_DOUBLE_EQ(1.0/log(2)+1/log(5),scores[i]);
        break;
      case 4:
        EXPECT_DOUBLE_EQ(1.0/log(2)+1/log(5),scores[i]);
        break;
      case 5:
        EXPECT_DOUBLE_EQ(1.0/log(5),scores[i]);
        break;
      case 6:
        EXPECT_DOUBLE_EQ(1.0/log(5),scores[i]);
        break;
      case 7:
        EXPECT_DOUBLE_EQ(1.0/log(5),scores[i]);
        break;
      case 8:
        EXPECT_DOUBLE_EQ(1.0/log(2),scores[i]);
        break;
    }
  }

  xfree(candidates);
  xfree(scores);
}

TEST_F(AdamicAdarTest, OneEdgeType) {
  for (int i=0; i < 100; i++) {
    int64_t timestamp = i+1;
    for (int j=i+1; j < 100; j++) {
      stinger_insert_edge_pair(S, 1, i, j, 1, timestamp);
    }
  }

  int64_t * candidates = NULL;
  double * scores = NULL;
  int64_t len = adamic_adar(S, 0, 0, &candidates, &scores);

  for (int64_t i = 0; i < len; i++) {
    switch (candidates[i]) {
      case 0:
        ADD_FAILURE() << "Start Vertex should not be a candidate";
        break;
      case 1:
        ADD_FAILURE() << "Vertex 1 should not be a candidate";
        break;
      case 2:
        ADD_FAILURE() << "Vertex 2 should not be a candidate";
        break;
      case 3:
        ADD_FAILURE() << "Vertex 3 should not be a candidate";
        break;
      case 4:
        EXPECT_DOUBLE_EQ(1.0/log(2)+1/log(5),scores[i]);
        break;
      case 5:
        EXPECT_DOUBLE_EQ(1.0/log(5),scores[i]);
        break;
      case 6:
        EXPECT_DOUBLE_EQ(1.0/log(5),scores[i]);
        break;
      case 7:
        EXPECT_DOUBLE_EQ(1.0/log(5),scores[i]);
        break;
      case 8:
        EXPECT_DOUBLE_EQ(1.0/log(2),scores[i]);
        break;
    }
  }

  xfree(candidates);
  xfree(scores);
}

int
main (int argc, char *argv[])
{
    printf("%d, %d, %lu \n", 1<<2, 1<<16, 1<<31);
    struct stinger_config_t * stinger_config;
    struct stinger * S;
    stinger_config = (struct stinger_config_t *)xcalloc(1,sizeof(struct stinger_config_t));
    stinger_config->nv = 1<<23;
    stinger_config->nebs = 1<<23;
    stinger_config->netypes = 2;
    stinger_config->nvtypes = 2;
    stinger_config->memory_size = 1<<45;
    S = stinger_new_full(stinger_config);
    xfree(stinger_config);

    printf("%d, %d, %lu \n", 1<<2, 1<<16, 1<<31);
    int64_t timestamp = 1;
	char *file_name = argv[1];
	if(!file_name) { fprintf(stderr, "File name could not be set.\n"); return 1; }
	FILE* f = fopen(argv[1], "rb");
	if(!f) { fprintf(stderr, "File %s could not be opened.\n", file_name); return 1; }
	int a, b;
	while(fscanf(f, " %i %i", &a, &b) == 2) {
		printf("A: %d B: %d\n", a, b);
      		stinger_insert_edge_pair(S, 0, a, b, 1, 1);
	}
	int64_t max_range = stinger_max_active_vertex(S);
	printf("max_range = %d\n", max_range);
	double start, end;
	double delta = 0;
	start = omp_get_wtime();
	printf("Start: %f \n", start);
        OMP("omp parallel for")
	for(int x = 0; x < max_range; x++)
	{
		int64_t *candidates = NULL;
		double *scores = NULL;
 		int64_t len = adamic_adar(S, x, 0, &candidates, &scores);
		free(candidates);
		free(scores);
	} 
	end = omp_get_wtime();
	printf("End: %f \n", end);
	delta = end - start;
	printf("Time: %f \n", delta);
}
