#include <stdio.h>
#include <math.h>

#define COMPUTE_EPSILON(type, eps_ptr) do {                \
    type eps = 1;                                          \
    int counter = 0;\
    while ((1 + eps / 2) != 1) {                           \
        eps /= 2;                                          \
        counter++;\
    }                                                      \
    *(eps_ptr) = eps;                                      \
} while (0)

#define COMPUTE_MANTISSA(type, mant_ptr) do {              \
    type one = 1;                                          \
    type val = 1;                                          \
    int bits = -1;                                         \
    while ((one + val) != one) {                           \
        val /= 2;                                          \
        bits++;                                            \
    }                                                      \
    *(mant_ptr) = bits;                                    \
} while (0)

#define COMPUTE_MIN_EXP(type, exp_ptr) do {                \
    type val = 1;                                          \
    int exp = 0;                                           \
    while (val / 2 != 0) {                                 \
        val /= 2;                                          \
        exp--;                                             \
    }                                                      \
    *(exp_ptr) = exp;                                      \
} while (0)

#define COMPUTE_MAX_EXP(type, exp_ptr, inf) do {           \
    type val = 1;                                          \
    int exp = 0;                                           \
    while (val * 2 != inf) {                               \
        val *= 2;                                          \
        exp++;                                             \
    }                                                      \
    *(exp_ptr) = exp;                                      \
} while (0)

#define COMPARE_VALUES(x, y) do { \
    printf("%s == %s → %s\n", #x, #y, ((x) == (y)) ? "true" : "false"); \
    printf("%s <  %s → %s\n", #x, #y, ((x) <  (y)) ? "true" : "false"); \
    printf("%s >  %s → %s\n", #x, #y, ((x) >  (y)) ? "true" : "false"); \
} while (0)

#define TEST_EPSILON(type) do {                       \
    type eps = 1;                                     \
    while ((1 + eps / 2) != 1) {                      \
        eps /= 2;                                     \
    }                                                 \
    type one = 1;                                     \
    type eps_half = eps / 2;                          \
    printf("Comparing %s values:\n", #type);          \
    COMPARE_VALUES(one, one + eps_half);              \
    COMPARE_VALUES(one, one + eps);                   \
    COMPARE_VALUES(one + eps, one + eps + eps_half);  \
} while (0)

int main() {
    float epsilon_float;
    double epsilon_double;

    COMPUTE_EPSILON(float, &epsilon_float);
    COMPUTE_EPSILON(double, &epsilon_double);

    int mantissa_float, mantissa_double;
    COMPUTE_MANTISSA(float, &mantissa_float);
    COMPUTE_MANTISSA(double, &mantissa_double);

    int min_exp_float, min_exp_double;
    COMPUTE_MIN_EXP(float, &min_exp_float);
    COMPUTE_MIN_EXP(double, &min_exp_double);

    int max_exp_float, max_exp_double;
    COMPUTE_MAX_EXP(float, &max_exp_float, INFINITY);
    COMPUTE_MAX_EXP(double, &max_exp_double, INFINITY);

    printf("=== float ===\n");
    printf("Машинное эпсилон: %.10e\n", epsilon_float);
    printf("Разрядность мантиссы: %d\n", mantissa_float);
    printf("Мин. экспонента: %d\n", min_exp_float);
    printf("Макс. экспонента: %d\n", max_exp_float);

    printf("\n=== double ===\n");
    printf("Машинное эпсилон: %.20e\n", epsilon_double);
    printf("Разрядность мантиссы: %d\n", mantissa_double);
    printf("Мин. экспонента: %d\n", min_exp_double);
    printf("Макс. экспонента: %d\n", max_exp_double);

    TEST_EPSILON(float);
    TEST_EPSILON(double);

    return 0;
}
