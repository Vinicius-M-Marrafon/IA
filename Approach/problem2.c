#include <stdio.h>

#define APPROACH_IMPLEMENTATION
#define SWAP_MUTATION

#include "approach.h"

/*
float graph[5][5] = {
    { 0,  5,  7,  4, 10},
    { 5,  0,  9,  3, 11},
    { 7,  9,  0, 13,  2},
    { 4,  3, 13,  0, 20},
    {10, 11,  2, 20,  0},
};*/

/*
[0.000000; 14.000000; 5.000000; 4.000000; 2.000000; 17.000000;
10.000000; 15.000000; 16.000000; 9.000000; 8.000000; 11.000000;
7.000000; 12.000000; 6.000000; 1.000000; 3.000000; 13.000000; 0.000000; ]
Fitness: 46.000000
*/
float graph[18][18] = {
    {  0,  5,  7,  4, 10,  1,  8,  6,  9, 12, 15, 11, 13,  3,  2, 14, 17, 16},
    {  5,  0,  9,  3, 11,  4,  2, 10,  7, 14, 13,  8, 16, 12,  6, 15, 17,  1},
    {  7,  9,  0, 13,  2,  5, 11,  8, 16,  3, 17, 12,  6, 14, 15, 10,  4,  1},
    {  4,  3, 13,  0, 20, 17, 12, 16, 15, 11,  8, 10,  5,  2,  7,  9,  1, 14},
    { 10, 11,  2, 20,  0,  1, 19, 18, 17, 16,  4,  6, 12,  9,  8, 15,  7,  5},
    {  1,  4,  5, 17,  1,  0,  7, 10,  3,  9, 12, 16, 11, 14,  2, 13,  8,  6},
    {  8,  2, 11, 12, 19,  7,  0,  9,  1, 13, 15,  5,  4, 16,  6,  3, 14, 10},
    {  6, 10,  8, 16, 18, 10,  9,  0, 14, 17,  2,  1,  3,  5, 15,  7, 12, 13},
    {  9,  7, 16, 15, 17,  3,  1, 14,  0,  5,  6,  2, 13, 18, 12, 11, 10,  4},
    { 12, 14,  3, 11, 16,  9, 13, 17,  5,  0, 19,  7,  8,  4, 10, 18,  2,  6},
    { 15, 13, 17,  8,  4, 12, 15, 11,  6, 19,  0, 18,  9,  1,  5,  3, 16,  2},
    { 11,  8, 12, 10,  6, 16,  5,  2,  2,  7, 18,  0, 14, 13,  4,  1,  9, 15},
    { 13, 16,  6,  5, 12, 11,  4,  3, 13,  8,  9, 14,  0, 20, 17, 10, 15,  7},
    {  3, 12, 14,  2,  9, 14, 16,  5, 18,  4,  1, 13, 20,  0,  7,  6, 11, 19},
    {  2,  6, 15,  7,  8,  2,  6, 15, 12, 10,  5,  4, 17,  7,  0,  9, 18, 13},
    { 14, 15, 10,  9, 15, 13,  3,  7, 11, 18,  3,  1, 10,  6,  9,  0,  5, 12},
    { 17, 17,  4,  1,  7,  8, 14, 12, 10,  2, 16,  9, 15, 11, 18,  5,  0,  3},
    { 16,  1,  1, 14,  5,  6, 10, 13,  4,  6,  2, 15,  7, 19, 13, 12,  3,  0}
};

float fitness(Gene *X)
{
    float sum = 0.0f;
    for (size_t i = 0; i < 18; i++) {
        size_t a, b;
        a = (size_t)((X + i)->code);
        b = (size_t)((X + i + 1)->code);
        sum += graph[a][b];
    }

    return sum;
}

int main(void)
{
    srand(time(0));

    Gene genes[] = {
        { .code = 0, .fixed = true  },
        { .code = 14, .fixed = false },
        { .code = 5, .fixed = false},
        { .code = 4, .fixed = false  },
        { .code = 2, .fixed = false },
        { .code = 17, .fixed = false },
        { .code = 10, .fixed = false },
        { .code = 15, .fixed = false },
        { .code = 16, .fixed = false },
        { .code = 9, .fixed = false },
        { .code = 8, .fixed = false },
        { .code = 11, .fixed = false },
        { .code = 7, .fixed = false },
        { .code = 12, .fixed = false },
        { .code = 6, .fixed = false },
        { .code = 1, .fixed = false },
        { .code = 3, .fixed = false },
        { .code = 13, .fixed = false },
        { .code = 0, .fixed = true },
    };   

    Ind model = {
        .count = 19,
        .high = 17,
        .low = 0,
        .genes = genes,
    };

    Spec s = spec_create(500, &model);
    float best_w = 1e3;
    Ind best = ind_alloc(19, 0, 17);

    for (size_t g = 0; g < 501; g++) {
        spec_select(&model, &s, fitness);
        float w = fitness(model.genes);
        
        if (w < best_w) {
            best_w = w;
            ind_clone(&best, &model);
        }

        spec_extermine(&s);
        s = spec_create(10, &best);
    }

    ind_print(&best);
    printf("Fitness: %f\n", best_w);
    ind_delete(&best);
    ind_delete(&model);
    spec_extermine(&s);
    
    return 0;
}