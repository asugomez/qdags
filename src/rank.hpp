#ifndef RANK_BV
#define RANK_BV

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include<bits/stdc++.h>

using namespace sdsl;
using namespace std;


class rank_bv_64 {
    uint64_t *seq;
    uint32_t *block;
    uint64_t u;  //bit vector length
    uint64_t n; // # ones

public:
    rank_bv_64() = default;

    rank_bv_64(bit_vector &bv) {
        uint64_t i;
        uint8_t byte_mask;
        uint32_t cur_word = 0, count = 0;

        u = bv.size();

        seq = new uint64_t[(u + 63) / 64]();
        block = new uint32_t[(u + 63) / 64]();

        for (i = 0; i < u; ++i) {

            if (i % 64 == 0)
                block[cur_word++] = count;

            if (bv[i]) { // count the 1s
                count++;
                seq[i / 64] |= (1L << (i % 64));
            } else
                seq[i / 64] &= ~(1L << (i % 64));

        }
        n = count;
    }

    // number of 1s in B[0,i]
    inline uint64_t rank(uint64_t i) {
        // 0x3f = 00111111
        // i >> 6 : dividir por 64 (tamaño bloque)
        // 0x3f : %64
        return block[i >> 6] + bits::cnt(seq[i >> 6] & ~(0xffffffffffffffff << (i & 0x3f)));
    }

    // los 4 bits que definen un nodo
    inline uint8_t get_4_bits(uint64_t start_pos) {
        return ((seq[start_pos >> 6] >> (start_pos & 0x3f)) & 0x0f);
    }

    inline uint8_t get_2_bits(uint64_t start_pos) {
        return ((seq[start_pos >> 6] >> (start_pos & 0x3f)) & 0x03);
    }

    /**
     * Count the number of 1s
     * @param start_pos
     * @return
     */
    inline uint64_t n_ones_4_bits(uint64_t start_pos) {
        // TODO: replace it by
        //uint64_t counter = bits::cnt(start_pos); // contar nro de reusltados
        uint8_t x = ((seq[start_pos >> 6] >> (start_pos & 0x3f)) & 0x0f);
        uint64_t counter = 0;
        for (int i = 0; i < 4; i++) {
            if (x & (1 << i)) {
                counter += 1;
            }
        }
        return counter;
    }

    inline uint64_t n_ones_2_bits(uint64_t start_pos) {
        // TODO: replace it by
        //uint64_t counter = bits::cnt(start_pos); // contar nro de reusltados
        uint8_t x = ((seq[start_pos >> 6] >> (start_pos & 0x3f)) & 0x03);
        uint64_t counter = 0;
        for (int i = 0; i < 2; i++) {
            if (x & (1 << i)) {
                counter += 1;
            }
        }
        return counter;
    }

    // number of bits in the bv
    inline uint64_t size() {
        return u;
    }

    inline uint64_t n_ones() {
        return n;
    }

    inline uint64_t size_in_bytes() {
        return sizeof(uint64_t) * ((u + 63) / 64) + sizeof(uint32_t) * (u + 63) / 64
               + sizeof(uint64_t *) + sizeof(uint32_t *)
               + 2 * sizeof(uint64_t);
    }

    void print_4_bits(uint64_t start_pos) {
        uint8_t x = ((seq[start_pos >> 6] >> (start_pos & 0x3f)) & 0x0f);

        for (int i = 0; i < 4; i++) {
            cout << ((x & (1 << i)) ? "1" : "0");
        }
        cout << " ";
    }

    void print_2_bits(uint64_t start_pos) {
        uint8_t x =  ((seq[start_pos >> 6] >> (start_pos & 0x3f)) & 0x03);

        for (int i = 0; i < 2; i++) {
            cout << ((x & (1 << i)) ? "1" : "0");
        }
        cout << " ";
    }

    void print(uint64_t k_d) {
        uint64_t i = 0;
        while (i < u) {
            if(k_d == 4){
                print_4_bits(i);
                i = i + 4;
            } else{
                print_2_bits(i);
                i = i + 2;
            }
        }
        cout << endl;
    }

};

#endif
