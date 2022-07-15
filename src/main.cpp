#include <random>
#include "cliqueEnum.h"
#include "polyEnum.h"

void random_graph(int n, double prob, int seed)
{
    assert(n > 0);
    assert(seed >= 0);
    int edges = 0;
    vector<vector<int>> adj;
    adj.resize(n);
    srand(seed);
    for (int i = 0; i < n; ++i) {
        for (int j = i+1;j < n; ++j) {
            if (rand() < RAND_MAX * prob) {
                adj[i].emplace_back(j);
                adj[j].emplace_back(i);
                edges += 2;
            }
        }
    }
    char outfile[1024];
    sprintf(outfile, "randomGraphs/random-%d-%f-%d.txt",n,prob,seed);
    FILE *out = fopen(outfile, "w");
    if (out == NULL) {
        printf("Failed to open %s\n", outfile);
        exit(1);
    }
    fprintf(out, "%d %d\n", n, edges/2);
    for (int i = 0; i < n; ++i) {
        for(auto v: adj[i]) {
            if (i < n) fprintf(out,"%d %d\n", i, v);
        }
    }
    fclose(out);

}

int main(int argc, char *argv[])
{
    Algorithm *mc = NULL;
    if (argc > 2 && strcmp(argv[1],"delay")==0) {
        mc = new PolyEnum();
    }
    else if (argc > 2 && strcmp(argv[1],"branch")==0) {
       mc = new CliqueEnum();
    }
    if (mc == NULL) return 1;
    int alg = 1;
    mc->read_graph(argv[2]);
    mc->setParameters(argc, argv);
    mc->run();
    delete mc;
    return 1;
}
