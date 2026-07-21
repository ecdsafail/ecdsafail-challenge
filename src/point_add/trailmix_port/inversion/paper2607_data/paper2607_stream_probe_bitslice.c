#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zstd.h>

enum {
    WIDTH = 581,
    WORK_WIDTH = 259,
    WORK1_START = 4,
    WORK2_START = 263,
    LT_START = 522,
    LQ_START = 531,
    LS_START = 540,
    LRP_START = 550,
    AUX_START = 559,
    LENGTH_WIDTH = 9,
    SHIFT_WIDTH = 10,
    SCHEDULE_STEPS = 1616,
    CASE_COUNT = 9,
    ALL_CASES = (1U << CASE_COUNT) - 1,
};

typedef uint16_t lane_t;

static const unsigned char P_LE[32] = {
    0x2f, 0xfc, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
};

static unsigned hex_nibble(char c);
static unsigned bit_p(unsigned i) { return (P_LE[i / 8] >> (i % 8)) & 1; }

static void parse_hex_le32(const char *text, unsigned char out[32]) {
    memset(out, 0, 32);
    const char *digits = text;
    while (*digits == '0' && digits[1] != '\0') digits++;
    size_t count = strlen(digits);
    assert(count > 0 && count <= 64);
    for (size_t digit = 0; digit < count; digit++) {
        unsigned value = hex_nibble(digits[count - 1 - digit]);
        unsigned byte = (unsigned)digit / 2;
        unsigned shift = 4 * ((unsigned)digit & 1);
        out[byte] |= (unsigned char)(value << shift);
    }
}

static int compare_le32(const unsigned char a[32], const unsigned char b[32]) {
    for (int i = 31; i >= 0; i--) {
        if (a[i] != b[i]) return a[i] > b[i] ? 1 : -1;
    }
    return 0;
}

static void subtract_le32(const unsigned char a[32], const unsigned char b[32],
                          unsigned char out[32]) {
    unsigned borrow = 0;
    for (unsigned i = 0; i < 32; i++) {
        unsigned sub = (unsigned)b[i] + borrow;
        out[i] = (unsigned char)((unsigned)a[i] - sub);
        borrow = (unsigned)a[i] < sub;
    }
    assert(borrow == 0);
}

static unsigned bit_length_le32(const unsigned char value[32]) {
    for (int byte = 31; byte >= 0; byte--) {
        if (value[byte] != 0) {
            unsigned top = value[byte];
            unsigned bits = 0;
            while (top != 0) {
                bits++;
                top >>= 1;
            }
            return 8U * (unsigned)byte + bits;
        }
    }
    return 0;
}

static unsigned lane_bit(const lane_t state[WIDTH], unsigned lane, unsigned case_index) {
    return (unsigned)((state[lane] >> case_index) & 1U);
}

static unsigned reg_value(const lane_t state[WIDTH], unsigned start, unsigned width,
                          unsigned case_index) {
    unsigned value = 0;
    for (unsigned i = 0; i < width; i++)
        value |= lane_bit(state, start + i, case_index) << i;
    return value;
}

static void apply_record(lane_t state[WIDTH], uint64_t word) {
    unsigned kind = word & 15;
    unsigned arity = (word >> 4) & 15;
    unsigned q0 = (word >> 8) & 1023;
    unsigned q1 = (word >> 18) & 1023;
    unsigned q2 = (word >> 28) & 1023;
    unsigned q3 = (word >> 38) & 1023;
    unsigned q4 = (word >> 48) & 1023;
    assert(q0 < WIDTH);
    if (kind == 1 && arity == 1) {
        state[q0] ^= ALL_CASES;
    } else if (kind == 2 && arity == 2) {
        assert(q1 < WIDTH);
        state[q1] ^= state[q0];
    } else if (kind == 3 && arity == 3) {
        assert(q1 < WIDTH && q2 < WIDTH);
        state[q2] ^= state[q0] & state[q1];
    } else if (kind == 7 && arity == 5) {
        assert(q1 < WIDTH && q2 < WIDTH && q3 < WIDTH && q4 < WIDTH);
        assert(state[q4] == 0);
        state[q3] ^= state[q0] & state[q1] & state[q2];
        assert(state[q4] == 0);
    } else {
        fprintf(stderr, "bad primitive %u/%u\n", kind, arity);
        exit(2);
    }
}

static unsigned char *decode(const char *path, size_t *size_out) {
    FILE *f = fopen(path, "rb");
    assert(f);
    fseek(f, 0, SEEK_END);
    long compressed_size = ftell(f);
    rewind(f);
    unsigned char *compressed = malloc((size_t)compressed_size);
    assert(compressed);
    assert(fread(compressed, 1, (size_t)compressed_size, f) == (size_t)compressed_size);
    fclose(f);
    ZSTD_DStream *stream = ZSTD_createDStream();
    assert(stream);
    size_t status = ZSTD_initDStream(stream);
    assert(!ZSTD_isError(status));
    size_t capacity = ZSTD_DStreamOutSize();
    unsigned char *decoded = malloc(capacity);
    assert(decoded);
    size_t used = 0;
    ZSTD_inBuffer input = { compressed, (size_t)compressed_size, 0 };
    while (input.pos < input.size) {
        if (capacity - used < ZSTD_DStreamOutSize()) {
            capacity *= 2;
            decoded = realloc(decoded, capacity);
            assert(decoded);
        }
        ZSTD_outBuffer output = { decoded + used, capacity - used, 0 };
        status = ZSTD_decompressStream(stream, &output, &input);
        assert(!ZSTD_isError(status));
        assert(output.pos != 0 || input.pos == input.size);
        used += output.pos;
    }
    assert(status == 0);
    ZSTD_freeDStream(stream);
    free(compressed);
    *size_out = used;
    return decoded;
}

static void run_chunks(lane_t state[WIDTH], int reverse, int argc, char **argv) {
    unsigned expected = reverse ? SCHEDULE_STEPS : 1;
    for (int k = 0; k < argc; k++) {
        int arg = reverse ? argc - 1 - k : k;
        size_t size;
        unsigned char *raw = decode(argv[arg], &size);
        assert(size >= 24 && !memcmp(raw, "P26EEA2\0", 8) && (size - 24) % 8 == 0);
        assert(*(uint32_t *)(raw + 8) == 256);
        assert(*(uint32_t *)(raw + 12) == WIDTH);
        unsigned start = *(uint32_t *)(raw + 16);
        unsigned end = *(uint32_t *)(raw + 20);
        if (reverse) {
            assert(end == expected);
            expected = start - 1;
        } else {
            assert(start == expected);
            expected = end + 1;
        }
        size_t count = (size - 24) / 8;
        for (size_t j = 0; j < count; j++) {
            size_t record = reverse ? count - 1 - j : j;
            uint64_t word;
            memcpy(&word, raw + 24 + 8 * record, 8);
            apply_record(state, word);
        }
        free(raw);
    }
    assert(expected == (reverse ? 0 : SCHEDULE_STEPS + 1));
}

static void rotate_right(lane_t *out, const lane_t *in, unsigned amount) {
    amount %= WORK_WIDTH;
    for (unsigned i = 0; i < WORK_WIDTH; i++) out[(i + amount) % WORK_WIDTH] = in[i];
}

static void print_work2_hex(const lane_t work2[WORK_WIDTH], unsigned case_index) {
    for (int nibble = 64; nibble >= 0; nibble--) {
        unsigned value = 0;
        for (unsigned bit = 0; bit < 4; bit++) {
            unsigned lane = (unsigned)nibble * 4 + bit;
            if (lane < WORK_WIDTH) value |= lane_bit(work2, lane, case_index) << bit;
        }
        printf("%x", value);
    }
}

static unsigned hex_nibble(char c) {
    if (c >= '0' && c <= '9') return (unsigned)(c - '0');
    if (c >= 'a' && c <= 'f') return (unsigned)(c - 'a' + 10);
    if (c >= 'A' && c <= 'F') return (unsigned)(c - 'A' + 10);
    fprintf(stderr, "invalid hex digit: %c\n", c);
    exit(2);
}

static void initialize_case(lane_t initial[WIDTH], unsigned case_index, const char *x_hex) {
    unsigned char x[32], half[32], used[32];
    lane_t mask = (lane_t)1U << case_index;
    parse_hex_le32(x_hex, x);
    unsigned carry = 0;
    for (int i = 31; i >= 0; i--) {
        unsigned combined = carry * 256U + P_LE[i];
        half[i] = (unsigned char)(combined >> 1);
        carry = combined & 1;
    }
    int high_half = compare_le32(x, half) > 0;
    if (high_half) {
        subtract_le32(P_LE, x, used);
        initial[2] |= mask;
    } else {
        memcpy(used, x, 32);
    }
    unsigned bitlen = bit_length_le32(used);
    assert(bitlen > 0);
    initial[WORK1_START] |= mask;
    for (unsigned bit = 0; bit < 256; bit++)
        if (bit_p(bit)) initial[WORK1_START + 258 - bit] |= mask;
    for (unsigned bit = 0; bit < 256; bit++)
        if ((used[bit / 8] >> (bit % 8)) & 1) initial[WORK2_START + 258 - bit] |= mask;
    for (unsigned i = 0; i < LENGTH_WIDTH; i++) initial[LQ_START + i] |= mask;
    for (unsigned i = 0; i < SHIFT_WIDTH; i++) initial[LS_START + i] |= mask;
    unsigned encoded = bitlen - 1;
    for (unsigned i = 0; i < LENGTH_WIDTH; i++)
        if ((encoded >> i) & 1) initial[LRP_START + i] |= mask;
}

static void inspect_case(const lane_t state[WIDTH], unsigned case_index, const char *x_hex) {
    assert(lane_bit(state, 0, case_index) == 0);
    assert(lane_bit(state, 1, case_index) == 0);
    assert(lane_bit(state, 3, case_index) == 0);
    for (unsigned i = 0; i < 22; i++)
        assert(lane_bit(state, AUX_START + i, case_index) == 0);
    assert(reg_value(state, LT_START, LENGTH_WIDTH, case_index) == 255);
    assert(reg_value(state, LQ_START, LENGTH_WIDTH, case_index) == 511);
    assert(reg_value(state, LRP_START, LENGTH_WIDTH, case_index) == 511);

    unsigned padding =
        (reg_value(state, LS_START, SHIFT_WIDTH, case_index) + 1) & 1023;
    lane_t canonical[WORK_WIDTH] = {0};
    rotate_right(canonical, state + WORK2_START, padding);
    printf("x=%s iter=%u padding=%u work2=", x_hex,
           lane_bit(state, 2, case_index), padding);
    print_work2_hex(canonical, case_index);
    printf("\n");
}

int main(int argc, char **argv) {
    assert(argc >= 2);
    const char *cases[CASE_COUNT] = {
        "1",
        "2",
        "3",
        "123456789abcdef",
        "6a09e667f3bcc908b2fb1366ea957d3e3adec17512775099da2f590a9c5d4a30",
        "5db3d742c265539d92ba16b83c5c1dc492ec1a6629ed23cc63905323d96efaef",
        "5db3d742c265539d92ba16b83c5c1dc492ec1a6629ed23cc63905323d950963b",
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2e",
        "d8f9f1d8b4f19c7a85c62a9b91f4eaa85283f694052036319488ed8fe28f2241",
    };
    lane_t initial[WIDTH] = {0};
    lane_t state[WIDTH];
    for (unsigned i = 0; i < CASE_COUNT; i++) initialize_case(initial, i, cases[i]);
    memcpy(state, initial, sizeof(state));
    run_chunks(state, 0, argc - 1, argv + 1);
    for (unsigned i = 0; i < CASE_COUNT; i++) inspect_case(state, i, cases[i]);
    run_chunks(state, 1, argc - 1, argv + 1);
    assert(!memcmp(state, initial, sizeof(state)));
    printf("PASS cases=%u forward_reverse=exact\n", CASE_COUNT);
    return 0;
}
