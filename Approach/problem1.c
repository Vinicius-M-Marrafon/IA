#define APPROACH_IMPLEMENTATION
#include "approach.h"


float fitness(Gene *X)
{
    float x = X->code;
    // A function with a loot of local max
    return sinf(x) + (2 * sinf(4 * x)) + (4 * sinf(8 * x));
}

int main(void) 
{
    srand(time(0));

    Ind model = ind_alloc(1, 0, 10);
    Spec s = spec_create(100, &model);
    s.sur_crit = MAXIMIZE;

    // History Best Fist
    float hbf = -1000.0f;

    // History Best Individual
    Ind hbi = ind_alloc(1, 0, 10);

    for (int g = 0; g < 20; g++) {
        Ind bgi = ind_alloc(1, 0, 10);
        
        spec_select(&bgi, &s, fitness);

        // Generation Best Fit
        float gbf = fitness(bgi.genes);
        if (gbf > hbf) {
            hbf = gbf;

            ind_clone(&hbi, &bgi);
        }

        printf("Generation: %d: Fitness = %f; Best Gene = [%f]\n", s.generation, gbf, bgi.genes[0].code);
        
        spec_extermine(&s);
        s = spec_create(100, &bgi);
        s.generation++;

        ind_delete(&bgi);
    }

    printf("=========================================\n");
    ind_print(&hbi);
    printf("fit = %f\n", hbf);
    
    spec_extermine(&s);
    ind_delete(&hbi);
    ind_delete(&model);

    return 0;
}