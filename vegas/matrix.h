//
//  matrix.h
//  vegas
//
//  Created by Lennart Oymanns on 19/06/16.
//  Copyright © 2016 Lennart Oymanns. All rights reserved.
//

#ifndef matrix_h
#define matrix_h
struct Matrix {
    double *m;
    int ncol;
    int nrow;
};

struct Matrix *matrix_new(int rows, int cols) {
    struct Matrix *m = malloc(sizeof(struct Matrix));
    m->ncol = cols;
    m->nrow = rows;
    m->m = malloc(sizeof(double) * rows * cols);
    int i;
    for (i = 0; i < rows * cols; i++) {
        m->m[i] = 0.0;
    }
    return m;
}

static void matrix_free(struct Matrix *m) {
    free(m->m);
    m->m = NULL;
    free(m);
    m = NULL;
}

static void matrix_set(struct Matrix *m, int row, int col, double v) {
    m->m[col + row * m->ncol] = v;
}

static double matrix_get(struct Matrix *m, int row, int col) {
    return m->m[col + row * m->ncol];
}

static void matrix_add_to_elem(struct Matrix *m, int row, int col, double v) {
    m->m[col + row * m->ncol] += v;
}

static void matrix_set_zero(struct Matrix *m) {
    int i;
    for (i = 0; i < m->nrow * m->ncol; i++) {
        m->m[i] = 0.0;
    }
}

static void matrix_add(struct Matrix *m, struct Matrix *p) {
    int i;
    for (i = 0; i < m->nrow * m->ncol; i++) {
        m->m[i] += p->m[i];
    }
}

static void matrix_mul_factor(struct Matrix *m, double v) {
    int i;
    for (i = 0; i < m->nrow * m->ncol; i++) {
        m->m[i] *= v;
    }
}

#endif /* matrix_h */
