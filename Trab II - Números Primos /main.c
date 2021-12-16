#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>    
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdlib.h>
#include <gmp.h>
#include <unistd.h>

/* 
 * Xorshift 128 bits. 
 * Taken from: https://en.wikipedia.org/wiki/Xorshift#Example_implementation
 */

typedef struct xorshift32_state {
  uint32_t value;
} xorshift32_state;

/* The state word must be initialized to non-zero */
uint32_t xorshift32(xorshift32_state *state) {
    /* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
    uint32_t x = state->value;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return state->value = x;
}

/* 
 * Park-Miller. 
 * Taken from: https://en.wikipedia.org/wiki/Xorshift#Example_implementation
 */

 uint32_t lcg_parkmiller(uint32_t *state) {
    uint64_t product = (uint64_t)*state * 48271;
    uint32_t x = (product & 0x7fffffff) + (product >> 31);

    x = (x & 0x7fffffff) + (x >> 31);
    return *state = x;
}

/* 
 * msleep(): Sleep for the requested number of milliseconds. 
 * Taken from https://stackoverflow.com/a/1157217/16808367
 */
int msleep(long msec) {
    struct timespec ts;
    int res;

    if (msec < 0) {
        errno = EINVAL;
        return -1;
    }

    ts.tv_sec = msec / 1e3;
    ts.tv_nsec = (msec % (long) 1e3) * 1e6;

    do {
        res = nanosleep(&ts, &ts);
    } while (res && errno == EINTR);

    return res;
}

/* 
 * Miller Rabin
 * Taken from https://www.sanfoundry.com/c-program-implement-rabin-miller-primality-test-check-number-prime/
 */ 

int miller_rabin_test(mpz_t p, int iterations) {
    mpz_t two;
    mpz_init_set_ui(two, 2);

    // Macro: int mpz_cmp_ui (const mpz_t op1, unsigned long int op2)
    // Compare op1 and op2. Return a positive value if op1 > op2, zero if op1 = op2, or a negative value if op1 < op2.
    if (mpz_cmp_ui(p, 2) < 0) return 0;

    mpz_t p_mod_2;
    mpz_init(p_mod_2);
    mpz_mod(p_mod_2, p, two);

    if (mpz_cmp_ui(p, 2) != 0 && mpz_cmp_ui(p_mod_2, 0) == 0) return 0;

    mpz_t y;
    mpz_init(y);
    mpz_sub_ui(y, p, 1);

    mpz_t y_mod_2;
    mpz_init(y_mod_2);

    // Divide y by 2 and assign the result to y only if y % 2 == 0
    while(1) {
        // Function: unsigned long int mpz_fdiv_r_ui (mpz_t r, const mpz_t n, unsigned long int d)
        mpz_fdiv_r_ui(y_mod_2, y, 2);

        if (mpz_cmp_ui(y_mod_2, 0) == 0) {
            mpz_fdiv_qr_ui(y, y_mod_2, y, 2);
        } else {
            break;
        }
    }

    // z = p - 1;
    mpz_t z;
    mpz_init(z);
    mpz_sub_ui(z, p, 1);

    for (int i = 0; i < iterations; i++) {
        mpz_t rand_num;
        mpz_init_set_ui(rand_num, rand());

        mpz_t a;
        mpz_init(a);

        // Function: void mpz_mod (mpz_t r, const mpz_t n, const mpz_t d)
        mpz_mod(a, rand_num, z);

        // Function: void mpz_add_ui (mpz_t rop, const mpz_t op1, unsigned long int op2)
        // Set rop to op1 + op2.
        mpz_add_ui(a, a, 1);

        mpz_t temp;
        mpz_init_set(temp, y);

        mpz_t mod;
        mpz_init(mod);

        // Function: void mpz_powm (mpz_t rop, const mpz_t base, const mpz_t exp, const mpz_t mod)
        // Set rop to (base raised to exp) modulo mod.
        mpz_powm(mod, a, temp, p);

        while (mpz_cmp(temp, z) != 0 && mpz_cmp_ui(mod, 1) != 0 && mpz_cmp(mod, z) != 0) {
            // Function: void mpz_mul (mpz_t rop, const mpz_t op1, const mpz_t op2)
            // Set rop to op1 times op2.
            mpz_mul(mod, mod, mod);
            mpz_mod(mod, mod, p);

            // Function: void mpz_mul_ui (mpz_t rop, const mpz_t op1, unsigned long int op2)
            // Set rop to op1 times op2.
            mpz_mul_ui(temp, temp, 2);
        }

        mpz_t temp_mod_2;
        mpz_init(temp_mod_2);
        mpz_mod(temp_mod_2, temp, two);

        if (mpz_cmp(mod, z) != 0 && mpz_cmp_ui(temp_mod_2, 0) == 0) {
            return 0;
        }
    }

    return 1;
}

/*
 * Fermat's test for checking primality
 */
int fermat_test(mpz_t p, int iterations) {
    // Macro: int mpz_cmp_ui (const mpz_t op1, unsigned long int op2)
    // Compare op1 and op2. Return a positive value if op1 > op2, zero if op1 = op2, or a negative value if op1 < op2.
    if (mpz_cmp_ui(p, 1) == 0) return 0;

    // z = p - 1;
    mpz_t z;
    mpz_init(z);
    mpz_sub_ui(z, p, 1);

    for (int i = 0; i < iterations; i++) {
        mpz_t rand_num;
        mpz_init_set_ui(rand_num, rand());

        mpz_t a;
        mpz_init(a);

        // Function: void mpz_mod (mpz_t r, const mpz_t n, const mpz_t d)
        mpz_mod(a, rand_num, z);

        // Function: void mpz_add_ui (mpz_t rop, const mpz_t op1, unsigned long int op2)
        // Set rop to op1 + op2.
        mpz_add_ui(a, a, 1);

        mpz_t mod;
        mpz_init(mod);

        // Function: void mpz_powm (mpz_t rop, const mpz_t base, const mpz_t exp, const mpz_t mod)
        // Set rop to (base raised to exp) modulo mod.
        mpz_powm(mod, a, z, p);

        if (mpz_cmp_ui(mod, 1) != 0) return 0;
    }

    return 1;
}

/* 
 * -------- MAIN ----------
 */

#define N_RAND_NUMBERS 11 // Generate numbers with 40, 56, 80, 128, 168, 224, 256, 512, 1024, 2048, 4096 bits.
#define UINT32_T_SIZE_IN_BITS 32 // Minimum of bits to generate a random number
#define BYTE_SIZE_IN_BITS 8 // Each position of a byte array
#define BYTES_IN_UINT32_T 4 // Number of bytes in uint32_t type

static const int n_bits[N_RAND_NUMBERS] = {40, 56, 80, 128, 168, 224, 256, 512, 1024, 2048, 4096};

// This is a random algorithm to generate a seed. Created by myself.
uint32_t generate_seed() {
    struct timeval first_tv, second_tv, third_tv;

    gettimeofday(&first_tv, NULL);
    const unsigned long long first_seed = (unsigned long long)(first_tv.tv_sec) * 1e3 + (unsigned long long)(first_tv.tv_usec) / 1e3; // millisecondsSinceEpoch

    gettimeofday(&second_tv, NULL);
    const unsigned long long second_seed = (unsigned long long)(second_tv.tv_sec) * 1e3 + (unsigned long long)(second_tv.tv_usec) / 1e3; // millisecondsSinceEpoch

    gettimeofday(&third_tv, NULL);
    const unsigned long long third_seed = (unsigned long long)(third_tv.tv_sec) * 1e3 + (unsigned long long)(third_tv.tv_usec) / 1e3; // millisecondsSinceEpoch

    return ((first_seed * second_seed) & third_seed) % UINT32_MAX;
}

void generate_rand(unsigned char *rand_bytes, int n_uint32_t, int n_chars, int generation_step, int is_xorshift) {
    // Generate each random 32 bits in a loop.
    for (int uint32_t_index = 0; uint32_t_index < n_uint32_t; uint32_t_index++) {
        xorshift32_state state;
        state.value = generate_seed();

        const uint32_t rand = is_xorshift > 0 ? xorshift32(&state) : lcg_parkmiller(&state.value);

        const unsigned char bytes[BYTES_IN_UINT32_T] = { (rand >> 24) & 0xFF, (rand >> 16) & 0xFF, (rand >> 8) & 0xFF, rand & 0xFF };
        const uint32_t temp_rand = (uint32_t) bytes[0] << 24 | (uint32_t) bytes[1] << 16 | (uint32_t) bytes[2] << 8 | (uint32_t) bytes[3];

        // Check if conversion from uint32_t to byte array and back worked.
        if (rand != temp_rand) {
            printf("[Error] Wrong conversion (Xorshift: %d) | uint32_t_index: %d | generation_step: %d | Expected: %lu | Got: %lu", is_xorshift, uint32_t_index, generation_step, (unsigned long) rand, (unsigned long) temp_rand);
            return exit(0);
        }

        // Make sure that it will have the numbers of bytes we want (defined by n_chars), since n_uint32_t can have more bytes than n_chars.
        const int offset = uint32_t_index * BYTES_IN_UINT32_T;
        const int size = offset + BYTES_IN_UINT32_T > n_chars ? n_chars - offset : BYTES_IN_UINT32_T;
        memcpy(&rand_bytes[offset], bytes, size * sizeof(unsigned char));
    }
}


void generate_random_number(char *rand_number, int generation_step, int is_xorshift, double* time_taken) {
    const int result_bits_size = n_bits[generation_step];
    int bits_size = result_bits_size;

    // Set "bits_size" to the closest number to "result_bits_size" that is greather than "result_bits_size" and divisible by 32.
    if (bits_size % UINT32_T_SIZE_IN_BITS != 0) {
        bits_size = UINT32_T_SIZE_IN_BITS;
        while (bits_size < result_bits_size) {
            bits_size += UINT32_T_SIZE_IN_BITS;
        }
    }

    const int n_uint32_t = bits_size / UINT32_T_SIZE_IN_BITS;
    const int n_chars = result_bits_size / BYTE_SIZE_IN_BITS;

    unsigned char rand_bytes[n_chars];

    // Measure how much time it took to generate each number;
    struct timeval start, end, result;

    gettimeofday(&start, NULL);
    generate_rand(rand_bytes, n_uint32_t, n_chars, generation_step, is_xorshift);
    gettimeofday(&end, NULL);   

    timersub(&end, &start, &result);

    *time_taken = (result.tv_sec * 1e6) + (result.tv_usec);

    // Converting unsigned char* to char*
    // Taken from https://stackoverflow.com/a/22260250/16808367O)
    for(int index = 0; index < n_chars; index++) {
        sprintf(rand_number + (index * 2), "%02x", rand_bytes[index]);
    }
    rand_number[n_chars * 2] = '\0';
}

int check_primality(mpz_t rand_num, double* time_taken, int is_miller_rabin) {
    // Measure how much time it took to check each number;
    struct timeval start, end, result;

    gettimeofday(&start, NULL);
    const int is_prime = is_miller_rabin > 0 ? miller_rabin_test(rand_num, 5) : fermat_test(rand_num, 10);
    gettimeofday(&end, NULL);

    timersub(&end, &start, &result);

    *time_taken = (result.tv_sec * 1e6) + (result.tv_usec);

    return is_prime;
}

int main(int argc, char **argv) {
    printf("\n---- STARTING ----\n\n");

    for (int step = 0; step < N_RAND_NUMBERS; step++) { 
        const int n_bits_step = n_bits[step];
        const int n_chars = n_bits_step / BYTE_SIZE_IN_BITS;

        printf("-------------- %d BITS --------------\n\n", n_bits_step);

        // rand_algorithm = 0 -> Parkmiller || rand_algorithm = 1 -> Xorshift
        for (int rand_algorithm = 0; rand_algorithm < 2; rand_algorithm++) {
            int rand_is_prime = 0;

            double rand_generation_time = 0;
            double miller_rabin_time = 0;
            double fermat_time = 0;
            double elapsed_time = 0;

            char rand_number[(2 * n_chars) + 1];

            mpz_t rand_num;
            mpz_init(rand_num);

            // Measure the number of iterations to find a rand prime number and the elapsed time (ms).
            int iterations = 0;
            struct timeval start, end, result;

            gettimeofday(&start, NULL);

            while (rand_is_prime < 1) {
                iterations++;
                if (iterations == INT_MAX - 2) iterations = 0;

                // printf("\rRunning iteration: %d | Max: %d | Until Max: %d", iterations, INT_MAX, INT_MAX - iterations);
                // fflush(stdout);

                // printf("\r");
                // fflush(stdout);

                // Generate a random number and copy it to rand_number byte array
                generate_random_number(rand_number, step, rand_algorithm, &rand_generation_time);

                // Use libgmp mpz_t type for big integers
                mpz_set_str(rand_num, rand_number, 16);

                // Run primality test for both algorithms, fermat first and then miller-rabin.
                rand_is_prime = check_primality(rand_num, &miller_rabin_time, 1);
                if (rand_is_prime > 0) {
                    rand_is_prime = check_primality(rand_num, &fermat_time, 0);
                }
            }

            gettimeofday(&end, NULL);
            timersub(&end, &start, &result);

            elapsed_time = ((result.tv_sec * 1e6) + (result.tv_usec)) / 1e3;

            const char *algorithm = rand_algorithm > 0 ? "Xorshift" : "Parkmiller";
            printf("[Found Prime] Time Taken(ms): %.2f || Iterations: %d\n\n", elapsed_time, iterations);
            printf("[%s] Random Number Generation Time Taken(us): %.2f\n", algorithm, rand_generation_time);
            printf("[Miller Rabin] Primality Check Time Taken(us): %.2f\n", miller_rabin_time);
            printf("[Fermat] Primality Check Time Taken(us): %.2f\n", fermat_time);
            gmp_printf("[%s] Rand Prime Number: %Zd\n\n", algorithm, rand_num);
            fflush(stdout);
        }
    }

    // check_primalitys(rand_numbers_xorshift, rand_numbers_parkmiller);
}

