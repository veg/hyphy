/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include "sequence.h"
#include "hy_strings.h"

_String nucl_alphabet = "ACGT-?",
        codon_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*?-.", full_alphabet,
        complete_nucl_alphabet = "AGCTUYRWSKMBDHVXN?0-.";

//______________________________________________________________________________
void initFullAlphabet(void) {
  _String f_a((unsigned long) 256);
  for (long i = 0; i < 256; i++) {
    f_a[i] = i;
  }
  full_alphabet = f_a;
}

//______________________________________________________________________________
unsigned int getNextCode(_String &s, long &p) {
  if (s.s_data[p] < 0) {
    unsigned int val = s.s_data[p++] & 0x7f;
    val *= 256;
    val += (unsigned char)(s.s_data[p++]);
    return val;
  }
  return s[p++];
}

//______________________________________________________________________________
unsigned char powers_of_2[9] = { 0, 2, 6, 14, 30, 62, 126, 254, 0 };
unsigned char real_powers_of_2[8] = { 1, 2, 4, 8, 16, 32, 64, 128 };

//______________________________________________________________________________
// auxiliary string bit write function
void writeBitsToString(_String &s, long &bit_at, char length_to_write) {
  long left_over = 8 - bit_at % 8, cur_pos = bit_at / 8;
  if (left_over >= length_to_write) { // will fit in current byte
    unsigned char value = (unsigned char) s[cur_pos];
    value += powers_of_2[left_over - 1] - powers_of_2[left_over - length_to_write];
    s[cur_pos] = value;
  } else {
    unsigned char value = (unsigned char) s[cur_pos];
    value += powers_of_2[left_over - 1] + 1;
    s[cur_pos] = value;
    char full_bytes = (length_to_write - left_over - 1) / 8;
    while (full_bytes) {
      s[++cur_pos] = 255;
      full_bytes--;
    }
    s[++cur_pos] = 254 - powers_of_2[8 - (length_to_write - left_over) % 8];
  }
  bit_at += length_to_write;
}
