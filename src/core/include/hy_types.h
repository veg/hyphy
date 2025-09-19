/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon42@uwo.ca)
  Steven Weaver (sweaver@ucsd.edu)

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

#ifndef __HY_TYPES__
#define __HY_TYPES__

#ifdef __ZLIB__
#include <zlib.h>
#endif

#include <stdio.h>

/**
 * @file hy_types.h
 * @brief Application-wide typedefs and enums.
 */

/** @brief Generic pointer type */
typedef void *hyPointer;
/** @brief Standard floating point type */
typedef double hyFloat;

/** @brief HyPhy Batch Language function types */
enum hyBLFunctionType {
  kBLFunctionAlwaysUpdate,
  kBLFunctionSkipUpdate,
  kBLFunctionLocal
};

/** @brief HyPhy Batch Language function argument types */
enum hyBLFunctionArgumentType {
  kBLFunctionArgumentNormal,
  kBLFunctionArgumentReference
};

/** @brief Tree definition phases */
enum hyTreeDefinitionPhase {
  kTreeNotBeingDefined,
  kTreeIsBeingParsed,
  kTreeNodeBeingCreated
};

/** @brief Comparison types */
enum hyComparisonType {
  kCompareLess = -1,
  kCompareEqual = 0,
  kCompareGreater = 1,
  kCompareUndefined = 0xff
};

/** @brief File open modes */
enum hyFileOpenMode {
  kFileRead,
  kFileReadBinary,
  kFileReadWrite,
  kFileWrite,
  kFileWriteBinary,
  kFileAppend
};

/** @brief A class for file operations */
class hyFile {
public:
  /**
   * @brief Construct a new hyFile object
   */
  hyFile(void) {
    _fileReference = NULL;
#ifdef __ZLIB__
    _fileReferenceDirect = NULL;
#endif
  }
  /**
   * @brief Open a file
   *
   * @param file_path The path to the file
   * @param mode The mode to open the file in
   * @param error Whether to throw an error if the file cannot be opened
   * @param compress Whether the file is compressed
   * @param buffer The buffer size
   * @return hyFile* A pointer to the opened file
   */
  static hyFile *openFile(const char *file_path, hyFileOpenMode mode,
                          bool error = false, bool compress = false,
                          long buffer = 1024 * 128);
  /**
   * @brief Lock the file
   */
  void lock(void);
  /**
   * @brief Unlock the file
   */
  void unlock(void);
  /**
   * @brief Rewind the file
   */
  void rewind(void);
  /**
   * @brief Seek to a position in the file
   *
   * @param offset The offset to seek to
   * @param whence The position to seek from
   */
  void seek(long, int);
  /**
   * @brief Flush the file
   */
  void flush(void);
  /**
   * @brief Close the file
   *
   * @return int The result of the close operation
   */
  int close();
  /**
   * @brief Check if the end of the file has been reached
   *
   * @return true if the end of the file has been reached, false otherwise
   */
  bool feof(void);
  /**
   * @brief Read from the file
   *
   * @param buffer The buffer to read into
   * @param size The size of each item to read
   * @param items The number of items to read
   * @return unsigned long The number of items read
   */
  unsigned long read(void *buffer, unsigned long size, unsigned long items);
  /**
   * @brief Write to the file
   *
   * @param buffer The buffer to write from
   * @param size The size of each item to write
   * @param count The number of items to write
   * @return size_t The number of items written
   */
  size_t fwrite(const void *buffer, size_t size, size_t count);
  /**
   * @brief Write a string to the file
   *
   * @param str The string to write
   * @return int The result of the write operation
   */
  int puts(const char *str);
  /**
   * @brief Write a character to the file
   *
   * @param chr The character to write
   * @return int The result of the write operation
   */
  int fputc(int chr);
  /**
   * @brief Get the current position in the file
   *
   * @return size_t The current position in the file
   */
  size_t tell(void);
  /**
   * @brief Read a character from the file
   *
   * @return int The character read from the file
   */
  int read_char(void);
#ifdef __ZLIB__
  /**
   * @brief Check if the file is valid
   *
   * @return true if the file is valid, false otherwise
   */
  inline bool valid(void) const {
    return _fileReference != NULL || _fileReferenceDirect != NULL;
  }
  gzFile _fileReference;
  FILE *_fileReferenceDirect;
  /**
   * @brief Check if the file is compressed
   *
   * @return true if the file is compressed, false otherwise
   */
  bool is_compressed(void) const { return _fileReference != NULL; }
#else
  /**
   * @brief Check if the file is valid
   *
   * @return true if the file is valid, false otherwise
   */
  inline bool valid(void) const { return _fileReference != NULL; }
  /**
   * @brief Check if the file is compressed
   *
   * @return true if the file is compressed, false otherwise
   */
  bool is_compressed(void) const { return false; }
  FILE *_fileReference;
#endif
};

#endif
