#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

/// TODO align everything to 100 columns
/// TODO add description to all parameters and return values of all procedures

/* Following implementation of MIDDLE is better than simple (left + right) / 2 for big values of
 * left and right, where sum of those is larger then MAX_INT.
 */
/// TODO use 1 >> instead of / 2 for division
#define MIDDLE(left, right) ((left) + ((right) - (left)) / 2)
/* Alternative */
//#define MIDDLE(left, right) (((left) + (right)) / 2)

/**
 * Binary search.
 *
 * @return if found, returns the position of found element
 *         if not found, returns the position of the smallest element bigger than key
 *
 * Complexity: O( log(N) )
 */
int
bs (const int *arr, int left, int right, int key, bool *found)
{
    int middle = MIDDLE(left, right);

    while (left <= right)
    {
        if (key < arr[middle])
            right = middle - 1;
        else if (key == arr[middle]) {
            *found = true;
            return middle;
        }
        else
            left = middle + 1;
        middle = MIDDLE(left, right);
    }

    *found = false;
    /* left points to the position of first bigger element */
    return left;
}

/**
 * Multi-key binary search - internal recursive call.
 */
static void
_mkbs (const int *arr, int arr_l, int arr_r, const int *keys, int M, int *results)
{
    /* end condition */
    if (M <= 0)
        return;

    int keys_middle = MIDDLE(0, M - 1);

    /* throw away half of keys, if the key on keys_middle is out */
    if (keys[keys_middle] < arr[arr_l]) {
        _mkbs(arr, arr_l, arr_r, &keys[keys_middle + 1], M - 1 - keys_middle, &results[keys_middle + 1]);
        return;
    }
    else if (keys[keys_middle] > arr[arr_r]) {
        _mkbs(arr, arr_l, arr_r, keys, keys_middle, results);
        return;
    }

    bool found;
    int pos = bs(arr, arr_l, arr_r, keys[keys_middle], &found);

    if (found)
        results[keys_middle] = pos;

    _mkbs(arr, arr_l, pos - 1, keys, keys_middle, results);
    _mkbs(arr, (found) ? pos + 1 : pos, arr_r, &keys[keys_middle + 1], M - 1 - keys_middle, &results[keys_middle + 1]);
}

/**
 * Multi-key binary search.
 *
 * Complexity (LaTex): O( \sum_{n = 0}^{log(M)}2^n*log(N/2^n) )
 */
void
mkbs (const int *arr, int N, const int *keys, int M, int *results)
{
    _mkbs(arr, 0, N - 1, keys, M, results);
}

struct boundaries {
    /* left boundary, points to the first element of an array where the result may occur */
    int l;
    /* right boundary, points to the last element of an array where the result may occur */
    int r;
};

/**
 * binary search wrapper of fmkbs, updates boundaries and writes result
 */
static void
_perform_bs (const int *arr,
             const int *keys,
             const int  M,
             const int  key_idx,
             struct boundaries *boundaries,
             int       *results)
{
    bool found;
    int pos = bs(arr, boundaries[key_idx].l, boundaries[key_idx].r, keys[key_idx], &found);  /// check ranges

    if (found)
        results[key_idx] = pos;
    if (key_idx > 0)
        boundaries[key_idx - 1].r = pos - 1;
    if (key_idx < M - 1)
        boundaries[key_idx + 1].l = (found) ? pos + 1 : pos;
}

// TODO selfsearch_performed dat na zaciatok, nech zabera menej miesta na stacku
/**
 * WIP: Experimental
 * TODO: reduce size & arguments
 *
 * Fast multi-key binary search - internal recursive call.
 */
static void
_fmkbs (const int *arr,
        const int  N,
        const int *keys,
        const int  M,
        const int  keys_l,
        const int  keys_r,
        const int  keys_middle,
        struct boundaries *boundaries,
        int       *results)
{
    const int m = keys_middle;
    const int arr_middle = MIDDLE(boundaries[m].l, boundaries[m].r);

    if (arr[arr_middle] == keys[m]) {    /* found */
        results[m] = arr_middle;

        if ((m - keys_l) > 0) {
            int l_half_m = MIDDLE(keys_l, m - 1);
            boundaries[l_half_m].r = (arr_middle > 0) ? arr_middle - 1 : 0;
            boundaries[l_half_m].l = boundaries[m].l;
            _fmkbs(arr, N, keys, M, keys_l, m - 1, l_half_m, boundaries, results);
        }
        if ((keys_r - m) > 0) {
            int r_half_m = MIDDLE(m + 1, keys_r);
            boundaries[r_half_m].l = (arr_middle < N - 1) ? arr_middle + 1 : N - 1;
            boundaries[r_half_m].r = boundaries[m].r;
            _fmkbs(arr, N, keys, M, m + 1, keys_r, r_half_m, boundaries, results);
        }

//        boundaries[m].l = boundaries[m].r = arr_middle; /// ??? TODO check if speeds up
    }
    else if (keys[m] < arr[arr_middle]) {    /* keys[m] will be somewhere on left side */
        int old_r = boundaries[m].r;
        bool selfsearch_performed = false;

        boundaries[m].r = (arr_middle > 0) ? arr_middle - 1 : arr_middle;

        /* at first search on the left side */
        if (!(keys[m] < arr[boundaries[m].l])) {
            if ((m - keys_l) > 0) {   // try to remove > 0
                int l_half_m = MIDDLE(keys_l, m - 1);
                boundaries[l_half_m].r = boundaries[m].r;
                boundaries[l_half_m].l = boundaries[m].l;
                _fmkbs(arr, N, keys, M, keys_l, m - 1, l_half_m, boundaries, results);
            }
            else {
                _perform_bs(arr, keys, M, m, boundaries, results);
                selfsearch_performed = true;
            }
        }

        /* then search on the right side */
        if ((keys_r - m) > 0) {
            int r_half_m = MIDDLE(m + 1, keys_r);
            boundaries[r_half_m].r = old_r;
            if (!selfsearch_performed)
                boundaries[r_half_m].l = boundaries[m].l;
            _fmkbs(arr, N, keys, M, m + 1, keys_r, r_half_m, boundaries, results);
        }
        else if (!selfsearch_performed) {
            _perform_bs(arr, keys, M, m, boundaries, results);
            selfsearch_performed = true;
        }

        if (!selfsearch_performed)
            _perform_bs(arr, keys, M, m, boundaries, results);
    }
    else {    /* keys[m] will be somewhere on right side */
        int old_l = boundaries[m].l;
        bool selfsearch_performed = false;

        boundaries[m].l = (arr_middle < N - 1) ? arr_middle + 1 : arr_middle;

        /* at first search on the right side */
        if (!(keys[m] > arr[boundaries[m].r])) {
            if ((keys_r - m) > 0) {
                int r_half_m = MIDDLE(m + 1, keys_r);
                boundaries[r_half_m].l = boundaries[m].l;
                boundaries[r_half_m].r = boundaries[m].r;
                _fmkbs(arr, N, keys, M, m + 1, keys_r, r_half_m, boundaries, results);
            }
            else {
                _perform_bs(arr, keys, M, m, boundaries, results);
                selfsearch_performed = true;
            }
        }

        /* then search on the left side */
        if ((m - keys_l) > 0) {
            int l_half_m = MIDDLE(keys_l, m - 1);
            boundaries[l_half_m].l = old_l;
            if (!selfsearch_performed)
                boundaries[l_half_m].r = boundaries[m].r;
            _fmkbs(arr, N, keys, M, keys_l, m - 1, l_half_m, boundaries, results);
        }
        else if (!selfsearch_performed) {
            _perform_bs(arr, keys, M, m, boundaries, results);
            selfsearch_performed = true;
        }

        if (!selfsearch_performed)
            _perform_bs(arr, keys, M, m, boundaries, results);
    }
}

/**
 * WIP: Experimental
 *
 * Fast multi-key binary search. (actually slower than mkbs)
 *
 * Complexity: O( M * (log((3 * N) / M) + log(M)) )   /// TODO check once more
 */
void
fmkbs (const int *arr, const int N, const int *keys, const int M, int *results)
{
    int middle_key = MIDDLE(0, M - 1);

    struct boundaries boundaries[M];
    /* set to defaults */
    boundaries[middle_key].l = 0;
    boundaries[middle_key].r = N - 1;

    _fmkbs(arr, N, keys, M, 0, M - 1, middle_key, boundaries, results);
}

/**
 * Bnary search called multiple times.
 *
 * Complexity: O( M * log(N) )
 */
void
m_times_bs (const int *arr, int N, const int *keys, int M, int *results)
{
    for (int i = 0; i < M; i++)
    {
        bool found;
        int pos = bs(arr, 0, N - 1, keys[i], &found);
        if (found)
            results[i] = pos;
    }
}

/***************************************************************************************************
 * GENERIC TESTSUITE
 */

#define TEST_OK    0
#define TEST_FAIL  1

/* qsort int comparison function */
int
int_cmp (const void *a, const void *b)
{
    /* casting pointer types */
    const int *ia = (const int *)a;
    const int *ib = (const int *)b;
    /* integer comparison: returns negative if b > a and positive if a > b */
    return *ia - *ib;
}

void
print_int_array (const int *a, int n)
{
    printf("[");
    for (int i = 0; i < n; i++)
        printf("%3d,", a[i]);
    printf(" ]\n");
}

/**
 * Check results correctness.
 *
 * @param[in] results_m_times_bs may be NULL, then check_results will compute it if needed
 */
int
check_results (const int *results_m_times_bs,
               const int *results_mkbs,
               const int *results_fmkbs,
               const int *arr,
               int N,
               const int *keys,
               int M,
               int test_num)
{
    if (!results_mkbs || !results_fmkbs) {
        printf(" No inputs: results_mkbs = %p, results_fmkbs = %p\n", results_mkbs, results_fmkbs);
        return TEST_FAIL;
    }
    else {
        /* compare mkbs & fmkbs results */
        if (0 != memcmp(results_mkbs, results_fmkbs, M * sizeof(int))) {
            /* results differ, compute the results also by m_times_bs */
            int *results_exp = NULL;
            if (!results_m_times_bs) {
                results_exp = (int *) malloc(M * sizeof(int));
                memset(results_exp, -1, M * sizeof(int));
                m_times_bs(arr, N, keys, M, results_exp);
            }

            printf(" Test %d FAILED --------------------------------------------\n", test_num);
            printf(" Expected results:\n");
            print_int_array((results_m_times_bs) ? results_m_times_bs : results_exp, M);
            printf(" fmkbs results:\n");
            print_int_array(results_fmkbs, M);
            printf(" mkbs results:\n");
            print_int_array(results_mkbs, M);
            printf(" arr[%d]:\n", N);
            print_int_array(arr, N);
            printf(" keys[%d]:\n", M);
            print_int_array(keys, M);
            putchar('\n');

            free(results_exp);
            return TEST_FAIL;
        }
    }

    return TEST_OK;
}

int
run_test (const int *arr, int N, const int *keys, int M, int test_num)
{
    int res;

    int *results_mkbs = (int *) malloc(M * sizeof(int));
    int *results_fmkbs = (int *) malloc(M * sizeof(int));

    memset(results_mkbs, -1, M * sizeof(int));
    memset(results_fmkbs, -1, M * sizeof(int));

    mkbs(arr, N, keys, M, results_mkbs);
    fmkbs(arr, N, keys, M, results_fmkbs);

    res = check_results(NULL, results_mkbs, results_fmkbs, arr, N, keys, M, test_num);

    free(results_mkbs);
    free(results_fmkbs);

    return res;
}

/* generic testsuite options */
#define ARR_VALUES_RANGE      20033 /* randomized search array values range */
#define KEY_VALUES_RANGE      27063 /* randomized key values range */
#define ARR_MAX_SIZE            999 /* (or N), it is randomized to values 0 .. ARR_MAX_SIZE */
#define KEYS_MAX_SIZE          1999 /* (or M), randomized too */
#define ALWAYS_GENERATE_DIFFERENT_TESTS


// TODO zbavit sa mallocov
void
run_generic_testsuite (const long test_cnt)
{
#ifdef ALWAYS_GENERATE_DIFFERENT_TESTS
    /* randomize seed */
    srand(time(NULL));
#endif

    int test_num = 0;
    int fail_cnt = 0;

    /* endless test */
    for (long a = 0; a < test_cnt; a++)
    {
        /** Create arr */
        /* with duplicates */
        int  N_dup        = (rand() % ARR_MAX_SIZE) + 1;
        int *arr_dup      = malloc(N_dup * sizeof(int));
        /* without duplicates */
        int  N            = 0;
        int *arr          = malloc(N_dup * sizeof(int));
        /* fill by random values */
        for (int i = 0; i < N_dup; i++)
            arr_dup[i] = rand() % ARR_VALUES_RANGE;
        qsort(arr_dup, N_dup, sizeof(int), int_cmp);
        /* remove duplicates */
        for (int i = 0; i < N_dup; i++)
            if (arr_dup[i] != arr_dup[i - 1] || i == 0)
                arr[N++] = arr_dup[i];
        free(arr_dup);

        /** Create keys */
        /* with duplicates */
        int  M_dup        = (rand() % KEYS_MAX_SIZE) + 1;
        int *keys_dup     = malloc(M_dup * sizeof(int));
        /* without duplicates */
        int  M            = 0;
        int *keys         = malloc(M_dup * sizeof(int));
        /* fill by random values */
        for (int i = 0; i < M_dup; i++)
            keys_dup[i] = rand() % KEY_VALUES_RANGE;
        qsort(keys_dup, M_dup, sizeof(int), int_cmp);
        /* remove duplicates */
        for (int i = 0; i < M_dup; i++)
            if (keys_dup[i] != keys_dup[i - 1] || i == 0)
                keys[M++] = keys_dup[i];
        free(keys_dup);

        fail_cnt += run_test(arr, N, keys, M, test_num);
        test_num++;

        free(arr);
        free(keys);
    }

    printf("tests ran: %ld, tests failed: %d\n", test_cnt, fail_cnt);
}

/***************************************************************************************************
 * MEASSUREMENTS
 */

#define CHECK_RESULTS_WHEN_MEASURING

#define t(func, ...) ({        \
    clock_t beg = clock();     \
    func(__VA_ARGS__);         \
    clock_t end = clock();     \
    ((long double) end - beg); \
})

/**
 * Do meassurements on pseudorandom arr and key values of size N, resp. M.
 */
void
meassure (const int N,
          const int M,
          const int arr_val_range,
          const int keys_val_range,
          long double *time_m_times_bs,
          long double *time_mkbs,
          long double *time_fmkbs)
{
    /** Create arr of size N */
    /* with duplicates */
    int *arr_dup      = malloc(2 * N * sizeof(int));
    /* without duplicates */
    int  arr_size     = 0;
    int *arr          = malloc(N * sizeof(int));
    while (arr_size < N)
    {
        arr_size = 0;
        /* fill by random values */
        for (int i = 0; i < 2 * N; i++)
            arr_dup[i] = rand() % arr_val_range;
        qsort(arr_dup, 2 * N, sizeof(int), int_cmp);
        /* remove duplicates */
        for (int i = 0; (i < 2 * N) && (arr_size < N); i++)
            if (arr_dup[i] != arr_dup[i - 1] || i == 0)
                arr[arr_size++] = arr_dup[i];
    }
    free(arr_dup);

    /** Create keys of size M */
    /* with duplicates */
    int *keys_dup     = malloc(2 * M * sizeof(int));
    /* without duplicates */
    int  keys_size    = 0;
    int *keys         = malloc(M * sizeof(int));
    while (keys_size < M)
    {
        keys_size = 0;
        /* fill by random values */
        for (int i = 0; i < 2 * M; i++)
            keys_dup[i] = rand() % keys_val_range;
        qsort(keys_dup, 2 * M, sizeof(int), int_cmp);
        /* remove duplicates */
        for (int i = 0; (i < 2 * M) && (keys_size < M); i++)
            if (keys_dup[i] != keys_dup[i - 1] || i == 0)
                keys[keys_size++] = keys_dup[i];
    }
    free(keys_dup);

    int *res_m_times_bs = malloc(M * sizeof(int));
    int *res_mkbs       = malloc(M * sizeof(int));
    int *res_fmkbs      = malloc(M * sizeof(int));

    memset(res_m_times_bs, -1, M * sizeof(int));
    memset(res_mkbs, -1, M * sizeof(int));
    memset(res_fmkbs, -1, M * sizeof(int));

    *time_m_times_bs = t(m_times_bs, arr, N, keys, M, res_m_times_bs);
    *time_mkbs       = t(mkbs,       arr, N, keys, M, res_mkbs);
    *time_fmkbs      = t(fmkbs,      arr, N, keys, M, res_fmkbs);

#ifdef CHECK_RESULTS_WHEN_MEASURING
    check_results(res_m_times_bs, res_mkbs, res_fmkbs, arr, N, keys, M, 0);
#endif

    free(res_m_times_bs);
    free(res_mkbs);
    free(res_fmkbs);
    free(arr);
    free(keys);
}

void
create_meassurements_csv (const char *csv_name)
{
    FILE *f = fopen(csv_name, "w");

    long double time_m_times_bs;
    long double time_mkbs;
    long double time_fmkbs;

    int meassurements_cnt = 20;

//    for (int M = 1000; M <= 400000; M += 1000)
    for (int N = 1000; N <= 400000; N += 1000)
    {
        long double sum_m_times_bs = 0;
        long double sum_mkbs = 0;
        long double sum_fmkbs = 0;

        /* meassure `meassurements_cnt` times for the same value of N resp. M and average it */
        for (int i = 0; i < meassurements_cnt; i++)
        {
//            meassure(200000, M, 50000000, 50000000, &time_m_times_bs, &time_mkbs, &time_fmkbs);
            meassure(N, 50000, 50000000, 50000000, &time_m_times_bs, &time_mkbs, &time_fmkbs);
            sum_m_times_bs += time_m_times_bs;
            sum_mkbs       += time_mkbs;
            sum_fmkbs      += time_fmkbs;
        }

        /* store average of all meassurements */
        fprintf(f, "%5d, %5d, %5d, %5d\n", N,
                                           (int) sum_m_times_bs / meassurements_cnt,
                                           (int) sum_mkbs       / meassurements_cnt,
                                           (int) sum_fmkbs      / meassurements_cnt);
    }

    fclose(f);
}

int
main ()
{
//    run_generic_testsuite(1000000);

    create_meassurements_csv("graph_N_var.csv");

    return 0;
}

