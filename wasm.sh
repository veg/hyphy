rm -rf CMakeFiles CMakeCache*

find tests/hbltests -name '*json' -delete
find tests/hbltests -name '*cache' -delete
find tests/hbltests -name '*fit' -delete
find tests/hbltests -name 'tempFile*' -delete
find tests/hbltests -name '*gard*' -delete
find tests/data ! -name '*.nex' -delete


emcmake cmake -DCMAKE_EXE_LINKER_FLAGS="-sTOTAL_STACK=2097152 -02 -sASSERTIONS=1 -sMODULARIZE=1 -sALLOW_MEMORY_GROWTH -sFORCE_FILESYSTEM=1 -sEXIT_RUNTIME=0 -s EXPORTED_RUNTIME_METHODS=["callMain","FS","PROXYFS","WORKERFS","UTF8ToString","getValue","AsciiToString"] -lworkerfs.js -lproxyfs.js -s INVOKE_RUN=0 -s ENVIRONMENT="web,worker" ${EM_FLAGS//-s /-s} -fwasm-exceptions --preload-file res@/hyphy --preload-file tests/hbltests@/tests"

   
emmake make -j hyphy

scp hyphy.data hyphy.js hyphy.wasm sergei@silverback.temple.edu:/data/shares/web/web/biowasm

