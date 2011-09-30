#include <stdio.h>
#include <string.h>

int main(int argc, char * argv[])
{
    FILE *fpi, *fpo;
    unsigned i = 0, j = 0;
    // is_line_start is true to start out with
    unsigned is_line_start = 1, is_slc = 0, is_mlc = 0;
    unsigned last_is_slash = 0, last_is_star = 0;
    char ibuf[1024], obuf[2048];

    if(argc != 3) {
        printf("usage: cl2hdr <CLFILE> <HEADERFILE>\n");
        return -1;
    }

    fpi = fopen(argv[1], "r");
    if(fpi < 0) {
        printf("Unable to open %s for reading.\n", argv[1]);
        return -1;
    }

    fpo = fopen(argv[2], "w");
    if(fpo < 0) {
        printf("Unable to open %s for writing.\n", argv[2]);
        return -1;
    }

    fprintf(fpo, "%s", "#define KERNEL_STRING \"");

    memset(&ibuf, 0, sizeof(ibuf));

    while(fread(&ibuf, sizeof(char), 1024, fpi) > 0) {
        j = 0;
        memset(&obuf, '\0', sizeof(obuf));
        for(i = 0; i < 1024 && ibuf[i] != '\0'; ++i) {
            // if we're a single-line comment
            if(last_is_slash && ibuf[i] == '/') {
                is_slc = 1;
                // only can be true if we're not a comment, which we are now
                last_is_slash = 0;
                continue;
            }

            // if we're a multi-line comment
            if(last_is_slash && ibuf[i] == '*') {
                is_mlc = 1;
                // only can be true if we're not a comment, which we are now
                last_is_slash = 0;
                continue;
            }

            if(last_is_slash && !(is_mlc || is_slc)) {
                obuf[j++] = '/';
            }

            /* last_is_star must also include is_mlc because
               it can only be turned on during is_mlc,
               and we don't print a '*' for last_is_star
               because of just this fact. */
            if(last_is_star && ibuf[i] == '/') {
                is_mlc = 0;
                continue;
            }

            // check this corner-case, if the last of the buf is a slash
            if(!(is_mlc || is_slc) && ibuf[i] == '/') {
                last_is_slash = 1;
                continue;
                // we'll print the slash if it's not a comment
            } else {
                last_is_slash = 0;
            }

            /* if we're a multi-line comment
               we need to check for the terminal part */
            if(is_mlc && ibuf[i] == '*') {
                last_is_star = 1;
                continue;
                // we'll print the star if it's not a comment
            } else {
                last_is_star = 0;
            }

            if(is_line_start && (ibuf[i] == '\t' || ibuf[i] == ' ')) {
                continue;
            }

            /* if it's a newline, escape it, add it, and continue
               ignoring comments, because we want to preserve line numbers */
            if (ibuf[i] == '\n') {
                // get rid of trailing whitespace
                for(; j > 1 && (obuf[j-1] == ' ' || obuf[j-1] == '\t'); --j);
                obuf[j++] = '\\';
                obuf[j++] = 'n';
                is_line_start = 1;
                is_slc = 0;
                continue;
            } else {
                is_line_start = 0;
            }

            // if it's a comment or we're at the start of a line
            // and we have whitespace
            if(is_slc || is_mlc) {
                continue;
            }

            // if it's a double-quote or backslash, escape it
            if(ibuf[i] == '"' || ibuf[i] == '\\') {
                obuf[j++] = '\\';
            }

            obuf[j++] = ibuf[i];
        }
        // make sure we terminate the output buffer properly
        obuf[j] = '\0';
        fprintf(fpo, "%s", obuf);
        memset(&ibuf, '\0', sizeof(ibuf));
    }

    fprintf(fpo, "%s", "\\n\"\n\n");

    return 0;
}
