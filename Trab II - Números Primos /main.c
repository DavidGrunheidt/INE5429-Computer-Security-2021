#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>    
#include <string.h>
#include <assert.h>
#include <stdlib.h>

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

    ts.tv_sec = msec / 1000;
    ts.tv_nsec = (msec % 1000) * 1000000;

    do {
        res = nanosleep(&ts, &ts);
    } while (res && errno == EINTR);

    return res;
}

/* 
 * -------- MAIN ----------
 */

#define N_RAND_NUMBERS 11 // Generate numbers with 40, 56, 80, 128, 168, 224, 256, 512, 1024, 2048, 4096 bits.
#define MAX_MS_TO_WAIT 10 // To generate different seeds
#define UINT32_T_SIZE_IN_BITS 32 // Minimum of bits to generate a random number
#define BYTE_SIZE_IN_BITS 8 // Each position of a byte array
#define BYTES_IN_UINT32_T 4 // Number of bytes in uint32_t type

static const int n_bits[N_RAND_NUMBERS] = {40, 56, 80, 128, 168, 224, 256, 512, 1024, 2048, 4096};

// This is a random algorithm to generate a seed. Created by myself.
uint32_t generate_seed() {
    struct timeval first_tv, second_tv, third_tv;

    gettimeofday(&first_tv, NULL);
    const unsigned long long first_seed = (unsigned long long)(first_tv.tv_sec) * 1000 + (unsigned long long)(first_tv.tv_usec) / 1000; // millisecondsSinceEpoch
    msleep(first_seed % MAX_MS_TO_WAIT);

    gettimeofday(&second_tv, NULL);
    const unsigned long long second_seed = (unsigned long long)(second_tv.tv_sec) * 1000 + (unsigned long long)(second_tv.tv_usec) / 1000; // millisecondsSinceEpoch
    msleep((first_seed * second_seed) % MAX_MS_TO_WAIT);

    gettimeofday(&third_tv, NULL);
    const unsigned long long third_seed = (unsigned long long)(third_tv.tv_sec) * 1000 + (unsigned long long)(third_tv.tv_usec) / 1000; // millisecondsSinceEpoch
    msleep(((first_seed * second_seed) & third_seed) % MAX_MS_TO_WAIT);

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

        const int offset = uint32_t_index * BYTES_IN_UINT32_T;
        memcpy(&rand_bytes[offset], bytes, BYTES_IN_UINT32_T * sizeof(char));
    }
}


void generate_random_number() {
    for (int generation_step = 0; generation_step < N_RAND_NUMBERS; generation_step++) { 
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
        const int n_chars = bits_size / BYTE_SIZE_IN_BITS;

        unsigned char rand_bytes_xorshift[n_chars];
        unsigned char rand_bytes_parkmiller[n_chars];

        // Measure how much time it took to generate each number;
        struct timeval start_xorshift, end_xorshift, start_parkmiller, end_parkmiller;

        gettimeofday(&start_xorshift, NULL);
        generate_rand(rand_bytes_xorshift, n_uint32_t, n_chars, generation_step, 1);
        gettimeofday(&end_xorshift, NULL);

        gettimeofday(&start_parkmiller, NULL);
        generate_rand(rand_bytes_parkmiller, n_uint32_t, n_chars, generation_step, 0);
        gettimeofday(&end_parkmiller, NULL);        

        int time_taken_xorshift = (end_xorshift.tv_sec + end_xorshift.tv_usec / 1e6 - start_xorshift.tv_sec - start_xorshift.tv_usec / 1e6) * 1000; // In milliseconds
        int time_taken_parkmiller = (end_parkmiller.tv_sec + end_parkmiller.tv_usec / 1e6 - start_parkmiller.tv_sec - start_parkmiller.tv_usec / 1e6) * 1000; // In milliseconds

        printf("[RGN %d] Needed: %d bits | Used: %d bits | n_uint32_t: %d | n_chars: %d\n", generation_step, result_bits_size, bits_size, n_uint32_t, n_chars);
        printf("[RGN %d] Time Taken(ms): %d | Rand Xorshift (hex): 0x", generation_step, time_taken_xorshift);
        for (int index = 0; index < n_chars; index++) {
            printf("%x", rand_bytes_xorshift[index]);
        }

        printf("\n[RGN %d] Time Taken(ms): %d | Rand Parkmiller (hex): 0x", generation_step, time_taken_parkmiller);
        for (int index = 0; index < n_chars; index++) {
            printf("%x", rand_bytes_parkmiller[index]);
        }
        printf("\n\n");
        fflush(stdout);
    }
}

int main(int argc, char **argv) {
    printf("\n---- STARTING ----\n\n");
    generate_random_number();
}

