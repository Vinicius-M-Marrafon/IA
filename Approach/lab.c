#include <stdio.h>

#define APPROACH_IMPLEMENTATION
#include "approach.h"

float fitness (Gene *X)
{
    return X->code;
}

int main(void)
{
    srand(time(0));
    Ind model = ind_alloc(1, -1, 1);
    Spec spec = spec_create(10, &model);
    
    spec_dump(&spec);

    Ind bgi = ind_alloc(1, -1, 1);
    spec_select(&bgi, &spec, fitness);

    ind_print(&bgi);

    spec.sur_crit = MAXIMIZE;
    spec_select(&bgi, &spec, fitness);

    ind_print(&bgi);

    spec_extermine(&spec);
    ind_delete(&model);
    ind_delete(&bgi);

    return 0;
}