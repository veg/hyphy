/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
   Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
   Art FY Poon    (apoon42@uwo.ca)
   Steven Weaver (sweaver@temple.edu)

Module Developers:
        Lance Hepler (nlhepler@gmail.com)
        Martin Smith (martin.audacis@gmail.com)

Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdf22@cam.ac.uk)

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

#pragma once
#include "avllist.h"
#include "avllistx.h"
#include "avllistxl.h"
#include "list.h"
#include "parser.h"
#include "simplelist.h"
#include "stdlib.h"

#define HY_TRANSLATION_TABLE_DNA 0x01
#define HY_TRANSLATION_TABLE_RNA 0x02
#define HY_TRANSLATION_TABLE_BINARY 0x04
#define HY_TRANSLATION_TABLE_PROTEIN 0x08


//_________________________________________________________
class _TranslationTable : public BaseObj {

private:
  static _List _list_of_default_tables;

public:
  _TranslationTable(void);
  _TranslationTable(unsigned char);
  _TranslationTable(_String &);
  /* 20100618: SLKP

          - new constructor (needed to handle ExecuteCase52 / Simulate properly)
            which takes an alphabet string and checks to see if it's a standard
     one DNA/RNA/Protein or Binary

   */
  _TranslationTable(_TranslationTable const &);
   
   _TranslationTable const& operator = (_TranslationTable const &);

  virtual ~_TranslationTable(void) {
    if (checkTable) {
      free(checkTable);
    }
  }
  virtual BaseRef makeDynamic(void) const;
  virtual void Duplicate(BaseRefConst);

  long TokenCode(char) const;
  char AmbigToLetter(long *, unsigned long) const;

  void AddBaseSet(_String const &);
  void SplitTokenCode(long, long *) const;

  long TokenResolutions(char token, long *buffer, bool gap_to_one = true) const;
  /**
    Given a character in this translation table, return the number of
    base alphabet characters that map to it, and populate buffer with
    their codes. For example `TokenResolutions ('T', buffer) will
    return 2 and set buffer[0] = 1, buffer[1] = 3, assuming the translation
    table has the standard IUPAC nucleotdie code

    @param token the character (unique or ambiguous) to translate
    @param buffer store the resolved characters (up to baseLength) here
    @param gap_to_one if `true`, map gaps (or equivalents) to fully ambiguous
    characters

    @return the number of resolutions

  */

  long MultiTokenResolutions(const _String &tokens, long *buffer,
                             bool gap_to_one = true) const;
  /**
   Given a string of several in this translation table, return the unique
   resolutions of characters to base, and populate buffer with their codes. For
   example `TokenResolutions ('ATR', buffer) will return 2 and set buffer[0] =
   12 [0*16+3*4+0], buffer[1] = 14 [0*16+3*4+2], assuming the translation table
   has the standard IUPAC nucleotdie code. Passing NULL as buffer will result in
   returning the code for the resolution (if UNIQUE), otherwise -1.

   @param tokens the characters (unique or ambiguous) to translate
   @param buffer store the resolved characters (up to baseLength) here [must be
   at least baseLength ^ length (token) long] can be set to NULL in which case
   the return behavior is modified
   @param gap_to_one if `true`, map gaps (or equivalents) to fully ambiguous
   characters

   @return the number of resolutions OR (if buffer == NULL) the code for the
   SINGLE resolution or -1 (multiple or invalid resolutions)
  */

  void AddTokenCode(char, _String const &);
  void PrepareForChecks(void);
  bool IsCharLegal(char);
  char GetSkipChar(void);
  char GetGapChar(void) const;
  const _String ConvertCodeToLetters(long, unsigned char) const;
  long LengthOfAlphabet(void) const;
  bool IsStandardBinary(void) const {
    return baseLength == 2 && baseSet.length() == 0;
  }
  bool IsStandardNucleotide(void) const {
    return baseLength == 4 && baseSet.length() == 0;
  }
  bool IsStandardAA(void) const {
    return baseLength == 20 && baseSet.length() == 0;
  }

  /** given a (possibly) ambiguous character
      expand it to a string to equivalent characters
      e.g. ExpandToken ('R') -> "AG" for IUPAC nucleotide data

      @param token the character to expand
      @return a string of complete expansions

   */
  const _String ExpandToken(char token) const;

  /** return the alphabet string for this table, so that code 0 maps to
   * result[0], etc
   *
   * @return the alphabet string (e.g. ACGT for standard DNA)
   */
  const _String &GetAlphabetString(void) const;


  _TranslationTable *MergeTables(_TranslationTable const *) const;
    
  bool  operator == (const _TranslationTable& ) const;

  static const _String &GetDefaultTable(long tableType);

  unsigned char baseLength;
  // number of "fundamental" tokens
  //(4 for nucl, ACGT; 20 for amino acids)

  _String tokensAdded, baseSet;

  _SimpleList translationsAdded;
  char *checkTable;
  // if null - then assume default translation table;
};

extern _TranslationTable hy_default_translation_table;
