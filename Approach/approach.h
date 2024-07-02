#ifndef APPROACH_H
#define APPROACH_H

#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
    MAXIMIZE,
    MINIMIZE
}
Crit;

typedef struct {
    float code;
    bool fixed;
}
Gene;

// Individual
typedef struct {
    Gene *genes;
    size_t count;
    float high, low;
}
Ind;

// Species
typedef struct {
    Ind *group;
    size_t count;
    size_t generation;
    // survival criterion
    Crit sur_crit;
}
Spec;

#define GENE_DESCRIBE(gene) gene_describe(#gene, gene)

float rand_float(void);

// ####################### GENE ###########################
Gene init_gene(float code, bool fixed);
void gene_mutate(Gene *dst, Gene *src, float eps);
void gene_print(Gene *gene);
void gene_describe(const char *label, Gene *gene);
// ####################### INDIVIDUAL ######################
Ind ind_alloc(size_t count, float low, float high);
void ind_mutate(Ind *dst, Ind *src, float eps);
void ind_clone(Ind *dst, Ind *src);
void ind_print(Ind *ind);
void ind_delete(Ind *ind);
// ####################### SPECIES #########################
Spec spec_create(size_t count, Ind *model);
void spec_select(Ind *dst, Spec *spec, float (*fitness)(Gene *));
void spec_extermine(Spec *spec);
void spec_dump(Spec *spec);

#ifdef APPROACH_IMPLEMENTATION

float rand_float()
{
    return ((float)rand()/ (float)RAND_MAX);
}

// ################################# GENE ########################################
Gene init_gene(float code, bool fixed)
{
    Gene gene;
    gene.code = code;
    gene.fixed = fixed;
    return gene;
}

void gene_mutate(Gene *dst, Gene *src, float eps)
{
    if (src != NULL && dst != NULL) {
        if (!dst->fixed)
            dst->code = src->code + eps;
    }
}

void gene_print(Gene *gene)
{
    printf("Gene: %f\n", gene->code);
}

void gene_describe(const char *label, Gene *gene)
{
    printf("%s = {\n", label);
    printf("\tFixed: %s\n", gene->fixed? "FIXED" : "CHANGEABLE");
    printf("\tGene: %f\n", gene->code);
    printf("}\n");
}

// ################################# INDIVIDUAL ########################################
Ind ind_alloc(size_t count, float low, float high)
{
    Ind ind;
    ind.count = count;
    ind.low = low;
    ind.high = high;
    ind.genes = malloc(count * sizeof(Gene));
    if (ind.genes != NULL) {
        for (size_t i = 0; i < count; i++) {
            float code = rand_float() * (high - low) + low;
            ind.genes[i] = init_gene(code, false);
        }

        return ind;
    }

    return (Ind){0};
}   

void ind_mutate(Ind *dst, Ind *src, float eps)
{
    if (dst != NULL && src != NULL) {
        if (dst->genes != NULL && src->genes != NULL && src->count == dst->count) {
            // Step mutation
#ifdef CODE_MUTATION
            for (size_t i = 0; i < dst->count; i++) {
                gene_mutate(&dst->genes[i], &src->genes[i], eps);
            }
#endif // CODE_MUTATION
#ifdef SWAP_MUTATION
            (void)eps;
            // Swap Mutation
            // The individual has at least two changeable genes?
            size_t changeable_genes = 0;
            for (size_t i = 0; i < src->count && changeable_genes < 2; i++) {
                if (!src->genes[i].fixed)
                    changeable_genes++;
            }

            if (changeable_genes >= 2) {
                size_t a, b;
                do 
                    a = rand() % src->count;
                while (src->genes[a].fixed);

                do {
                    b += (a + (size_t)((rand_float()) * (2 * src->count) - src->count));
                    b %= src->count;
                }
                while (src->genes[b].fixed || a == b);

                // Swap;
                // printf("Swapping %u x %u\n", a, b);
                Gene c = dst->genes[a];
                dst->genes[a] = src->genes[b];
                dst->genes[b] = c;
            }

#endif // SWAP_MUTATION
        }
    }
}

void ind_clone(Ind *dst, Ind *src)
{
    if (dst != NULL && src != NULL) {
        if (dst->genes != NULL && src->genes != NULL && dst->count == src->count) {
            for (size_t i = 0; i < dst->count; i++) {
                dst->genes[i] = src->genes[i];
            }
        }
    }
}

void ind_print(Ind *ind)
{
    if (ind != NULL) {
        printf("[");
        for (size_t i = 0; i < ind->count; i++) {
            printf("%f; ", ind->genes[i].code);
        }
        printf("]\n");
    }
}

void ind_delete(Ind *ind)
{
    if (ind != NULL) { 
        if (ind->genes != NULL) { 
            free(ind->genes);
            ind->count = 0;
        }
    }
}

// ################################# SPECIES ########################################
Spec spec_create(size_t count, Ind *model)
{
    Spec spec;
    spec.count = count;
    spec.generation = 0;
    spec.sur_crit = MINIMIZE;
    spec.group = malloc(count * sizeof(Ind));

    if (spec.group != NULL) {
        for (size_t i = 0; i < count; i++) {
            spec.group[i] = ind_alloc(model->count, model->low, model->high);
            ind_clone(&spec.group[i], model);
            if (i > 0) 
                ind_mutate(&spec.group[i], model, rand_float());
        }

        return spec;
    }

    return (Spec){0};
}

// Select the best individual of the specie using fitness function
void spec_select(Ind *dst, Spec *spec, float (*fitness)(Gene *))
{
    if (dst != NULL && spec != NULL && fitness != NULL) {
        if (dst->genes != NULL && spec->group != NULL) {
            float best_fitness;
            size_t best_index = 0;
            
            if (spec->sur_crit == MINIMIZE) best_fitness = 1e5;
            else if (spec->sur_crit == MAXIMIZE) best_fitness = -1e5;

            for (size_t i = 0; i < spec->count; i++) {
                float ind_fitness = fitness(spec->group[i].genes);
                
                // MINIMIZE
                if (spec->sur_crit == MINIMIZE) {
                    if (ind_fitness < best_fitness) {
                        best_fitness = ind_fitness;
                        best_index = i; 
                    }
                }

                // MAXIMIZE
                else {
                    if (ind_fitness > best_fitness) {
                        best_fitness = ind_fitness;
                        best_index = i; 
                    }
                }
            }

            ind_clone(dst, &spec->group[best_index]);
        }
    }
}

void spec_extermine(Spec *spec)
{
    if (spec != NULL) {
        if (spec->group != NULL) {
            for (size_t i = 0; i < spec->count; i++) {
                ind_delete(&spec->group[i]);
            }

            free(spec->group);
            spec->group = NULL;
            spec->count = 0;
        }
    }
}

void spec_dump(Spec *spec) {
    if (spec != NULL) {
        if (spec->group != NULL) {
            printf("SPECIE: {\n");
            for (size_t i = 0; i < spec->count; i++) {
                printf("\t");
                ind_print(&spec->group[i]);
            }
            printf("} END SPECIE\n");
        }
    }
}

#endif // APPRAOCH_IMPLEMENTATION

#endif // APPROACH_H