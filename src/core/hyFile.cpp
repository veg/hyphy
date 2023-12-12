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

#include "hy_types.h"
#include "hy_strings.h"
#include "global_things.h"

//____________________________________________________________________________________
    hyFile* hyFile::openFile (const char * file_path, hyFileOpenMode mode , bool error, bool compress, long buffer) {
        
        _String str_mode;
        auto handle_error = [&](void * ptr)->void {
            if (!ptr && error) {
                hy_global::HandleApplicationError (_String("Could not open file '") & *file_path & "' with mode '" & str_mode & "'.");
            }
        };
        
        hyFile* f = new hyFile;
#ifdef __ZLIB__
        if (file_path) {
            
            _String str_mode;
            switch (mode) {
                case kFileRead:
                    str_mode = "r";
                    break;
                case kFileReadBinary:
                    str_mode = "rb";
                    break;
                case kFileWrite:
                    str_mode = compress ? "w" : "wT";
                    break;
                case kFileWriteBinary:
                    str_mode = compress ? "wb" : "wbT";
                    break;
               case kFileAppend:
                    str_mode = compress ? "a" : "aT";
                    break;
                case kFileReadWrite:
                    str_mode = "w+";
                    break;
            }
            
            if (mode != kFileReadWrite) {
                f->_fileReference = gzopen (file_path, (char const*)str_mode);
                if (f->_fileReference) {
                    if (gzdirect(f->_fileReference)) {
                        gzclose (f->_fileReference);
                        f->_fileReference = NULL;
                        f->_fileReferenceDirect = ::fopen (file_path, (const char*)str_mode);
                        handle_error (f->_fileReferenceDirect);
                    }
                } else {
                    handle_error (f->_fileReference);
                }
            } else {
                f->_fileReferenceDirect = ::fopen (file_path, (const char*)str_mode);
                handle_error (f->_fileReferenceDirect);
            }
            if (!f->valid()) {
                delete f;
                f = nil;
            }
}
#else
        if (file_path) {
            _String str_mode;
            switch (mode) {
                case kFileRead:
                    str_mode = "r";
                    break;
                case kFileReadBinary:
                    str_mode = "rb";
                    break;
                case kFileWrite:
                    str_mode = "w";
                    break;
                case kFileWriteBinary:
                    str_mode = "wb";
                    break;
                case kFileAppend:
                    str_mode = "a";
                    break;
                case kFileReadWrite:
                    str_mode = "w+";
                    break;
            }
            
            f->_fileReference = ::fopen (file_path, (const char*)str_mode);
            handle_error (f->_fileReference);
        }
        if (!f->_fileReference ) {
            delete f;
            f = nil;
        }
#endif
        
        return f;
        
    }

    //____________________________________________________________________________________
    int hyFile::close (void) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    return fclose (_fileReferenceDirect);
                } else {
                    return gzclose (_fileReference);
                }
            #else
                return fclose (_fileReference);
            #endif
            _fileReference = NULL;
        }
        return 0;
    }

    //____________________________________________________________________________________
    void hyFile::lock (void) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    flockfile (_fileReferenceDirect);
                }
                //gzclose (_fileReference);
            #else
                flockfile (_fileReference);
            #endif
        }
    }

    //____________________________________________________________________________________
    void hyFile::flush (void) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    fflush (_fileReferenceDirect);
                } else {
                    gzflush (_fileReference, Z_SYNC_FLUSH);
                }
            #else
                fflush(_fileReference);
            #endif
        }
    }

    //____________________________________________________________________________________
    void hyFile::unlock (void) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    funlockfile (_fileReferenceDirect);
                }
            #else
                funlockfile (_fileReference);
            #endif
        }
    }

    //____________________________________________________________________________________
    void hyFile::rewind (void) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    ::rewind (_fileReferenceDirect);
                } else {
                    gzrewind (_fileReference);
                }
            #else
                ::rewind (_fileReference);
            #endif
        }
    }

    //____________________________________________________________________________________
    void hyFile::seek (long pos, int whence) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    fseek (_fileReferenceDirect, pos, whence);
                } else {
                    gzseek (_fileReference, pos, whence);
                }
            #else
                fseek (_fileReference, pos, whence);
            #endif
        }
    }
    
    //____________________________________________________________________________________

    size_t hyFile::tell (void) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    return ftell (_fileReferenceDirect);
                } else {
                    return gztell (_fileReference);
                }
            #else
                return ftell (_fileReference);
            #endif
        }
        return 0;
    }
    //____________________________________________________________________________________
    bool hyFile::feof (void) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    return feof_unlocked (_fileReferenceDirect);
                } else {
                    return gzeof (_fileReference);
                }
            #else
                return feof_unlocked (_fileReference);
            #endif
        }
        return true;
    }
    //____________________________________________________________________________________
    int hyFile::puts( const char* buffer) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    return ::fputs (buffer, _fileReferenceDirect);
                } else {
                    return gzputs (_fileReference, buffer);
                }
           
            #else
                return ::fputs (buffer, _fileReference);
            #endif
        }
        return -1;
    }

    //____________________________________________________________________________________
    int hyFile::fputc( int chr) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    return ::fputc (chr, _fileReferenceDirect);
                } else {
                    return gzputc (_fileReference, chr);
                }
            #else
                return ::fputc (chr, _fileReference);
            #endif
        }
        return -1;
    }

    //____________________________________________________________________________________
    size_t hyFile::fwrite( const void* buffer, size_t size, size_t count) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    return ::fwrite (buffer, size, count, _fileReferenceDirect);
                } else {
                    return gzfwrite (buffer, size, count, _fileReference);
                }
            #else
                return ::fwrite (buffer, size, count, _fileReference);
            #endif
        }
        return -1;
    }
    //____________________________________________________________________________________
    int hyFile::getc (void) {
        if (valid()) {
            #ifdef __ZLIB__
                if (_fileReferenceDirect) {
                    return getc_unlocked (_fileReferenceDirect);
                } else {
                    return gzgetc (_fileReference);
                }
               
            #else
                return getc_unlocked (_fileReference);
            #endif
        }
        return 0;
    }

    //____________________________________________________________________________________
    unsigned long hyFile::read (void* buffer, unsigned long size, unsigned long items) {
        if (valid()) {
        #ifdef __ZLIB__
            if (_fileReferenceDirect) {
                return ::fread (buffer, size, items, _fileReferenceDirect);
            } else {
                return gzfread (buffer, size, items, _fileReference);
            }
            
        #else
            return ::fread (buffer, size, items, _fileReference);
        #endif
        }
        return 0;
    }
