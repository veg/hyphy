// This code implements the `-sMODULARIZE` settings by taking the generated
// JS program code (INNER_JS_CODE) and wrapping it in a factory function.

// Single threaded MINIMAL_RUNTIME programs do not need access to
// document.currentScript, so a simple export declaration is enough.
var Module = (() => {
  // When MODULARIZE this JS may be executed later,
  // after document.currentScript is gone, so we save it.
  // In EXPORT_ES6 mode we can just use 'import.meta.url'.
  var _scriptName =
      typeof document != 'undefined' ? document.currentScript?.src : undefined;
  return async function(moduleArg = {}) {
    var moduleRtn;

    // include: shell.js
    // The Module object: Our interface to the outside world. We import
    // and export values on it. There are various ways Module can be used:
    // 1. Not defined. We create it here
    // 2. A function parameter, function(moduleArg) => Promise<Module>
    // 3. pre-run appended it, var Module = {}; ..generated code..
    // 4. External script tag defines var Module.
    // We need to check if Module already exists (e.g. case 3 above).
    // Substitution will be replaced with actual code on later stage of the
    // build, this way Closure Compiler will not mangle it (e.g. case 4. above).
    // Note that if you want to run closure, and also to use Module
    // after the generated code, you will need to define   var Module = {};
    // before the code. Then that object will be used in the code, and you
    // can continue to use Module afterwards as well.
    var Module = moduleArg;

    // Determine the runtime environment we are in. You can customize this by
    // setting the ENVIRONMENT setting at compile time (see settings.js).

    // Attempt to auto-detect the environment
    var ENVIRONMENT_IS_WEB = typeof window == 'object';
    var ENVIRONMENT_IS_WORKER = typeof WorkerGlobalScope != 'undefined';
    // N.b. Electron.js environment is simultaneously a NODE-environment, but
    // also a web environment.
    var ENVIRONMENT_IS_NODE = typeof process == 'object' &&
                              process.versions?.node &&
                              process.type != 'renderer';
    var ENVIRONMENT_IS_SHELL =
        !ENVIRONMENT_IS_WEB && !ENVIRONMENT_IS_NODE && !ENVIRONMENT_IS_WORKER;

    // --pre-jses are emitted after the Module integration code, so that they
    // can refer to Module (if they choose; they can also define Module)
    // include: /var/folders/_6/5pbbnj3j5x1b_vnpy1q7fzlc0000gp/T/tmphcfvyuxw.js

    Module['expectedDataFileDownloads'] ??= 0;
    Module['expectedDataFileDownloads']++;
    (() => {
      // Do not attempt to redownload the virtual filesystem data when in a
      // pthread or a Wasm Worker context.
      var isPthread = typeof ENVIRONMENT_IS_PTHREAD != 'undefined' &&
                      ENVIRONMENT_IS_PTHREAD;
      var isWasmWorker = typeof ENVIRONMENT_IS_WASM_WORKER != 'undefined' &&
                         ENVIRONMENT_IS_WASM_WORKER;
      if (isPthread || isWasmWorker)
        return;
      async function loadPackage(metadata) {

        var PACKAGE_PATH = '';
        if (typeof window === 'object') {
          PACKAGE_PATH = window['encodeURIComponent'](
              window.location.pathname.substring(
                  0, window.location.pathname.lastIndexOf('/')) +
              '/');
        } else if (typeof process === 'undefined' &&
                   typeof location !== 'undefined') {
          // web worker
          PACKAGE_PATH =
              encodeURIComponent(location.pathname.substring(
                                     0, location.pathname.lastIndexOf('/')) +
                                 '/');
        }
        var PACKAGE_NAME = 'hyphy.data';
        var REMOTE_PACKAGE_BASE = 'hyphy.data';
        var REMOTE_PACKAGE_NAME =
            Module['locateFile']?.(REMOTE_PACKAGE_BASE, '') ??
            REMOTE_PACKAGE_BASE;
        var REMOTE_PACKAGE_SIZE = metadata['remote_package_size'];

        async function fetchRemotePackage(packageName, packageSize) {

          Module['dataFileDownloads'] ??= {};
          try {
            var response = await fetch(packageName);
          } catch (e) {
            throw new Error(`Network Error: ${packageName}`, {e});
          }
          if (!response.ok) {
            throw new Error(`${response.status}: ${response.url}`);
          }

          const chunks = [];
          const headers = response.headers;
          const total = Number(headers.get('Content-Length') ?? packageSize);
          let loaded = 0;

          Module['setStatus']?.('Downloading data...');
          const reader = response.body.getReader();

          while (1) {
            var {done, value} = await reader.read();
            if (done)
              break;
            chunks.push(value);
            loaded += value.length;
            Module['dataFileDownloads'][packageName] = {loaded, total};

            let totalLoaded = 0;
            let totalSize = 0;

            for (const download of Object.values(Module['dataFileDownloads'])) {
              totalLoaded += download.loaded;
              totalSize += download.total;
            }

            Module['setStatus']?.(
                `Downloading data... (${totalLoaded}/${totalSize})`);
          }

          const packageData = new Uint8Array(
              chunks.map((c) => c.length).reduce((a, b) => a + b, 0));
          let offset = 0;
          for (const chunk of chunks) {
            packageData.set(chunk, offset);
            offset += chunk.length;
          }
          return packageData.buffer;
        }

        var fetchedCallback;
        var fetched = Module['getPreloadedPackage']?.(REMOTE_PACKAGE_NAME,
                                                      REMOTE_PACKAGE_SIZE);

        if (!fetched) {
          // Note that we don't use await here because we want to execute the
          // the rest of this function immediately.
          fetchRemotePackage(REMOTE_PACKAGE_NAME, REMOTE_PACKAGE_SIZE)
              .then((data) => {
                if (fetchedCallback) {
                  fetchedCallback(data);
                  fetchedCallback = null;
                } else {
                  fetched = data;
                }
              });
        }

        async function runWithFS(Module) {

          function assert(check, msg) {
            if (!check)
              throw new Error(msg);
          }
          Module['FS_createPath']("/", "hyphy", true, true);
          Module['FS_createPath']("/hyphy", "GeneticCodes", true, true);
          Module['FS_createPath']("/hyphy", "SubstitutionClasses", true, true);
          Module['FS_createPath']("/hyphy/SubstitutionClasses", "AAEFV", true,
                                  true);
          Module['FS_createPath']("/hyphy/SubstitutionClasses", "CodonEFV",
                                  true, true);
          Module['FS_createPath']("/hyphy/SubstitutionClasses", "Heterogeneity",
                                  true, true);
          Module['FS_createPath']("/hyphy/SubstitutionClasses", "NucEFV", true,
                                  true);
          Module['FS_createPath']("/hyphy", "SubstitutionModels", true, true);
          Module['FS_createPath']("/hyphy/SubstitutionModels", "Aminoacid",
                                  true, true);
          Module['FS_createPath']("/hyphy/SubstitutionModels", "Binary", true,
                                  true);
          Module['FS_createPath']("/hyphy/SubstitutionModels", "Codon", true,
                                  true);
          Module['FS_createPath']("/hyphy/SubstitutionModels", "Nucleotide",
                                  true, true);
          Module['FS_createPath']("/hyphy/SubstitutionModels", "User", true,
                                  true);
          Module['FS_createPath']("/hyphy/SubstitutionModels/User",
                                  "Nucleotide", true, true);
          Module['FS_createPath']("/hyphy", "TemplateBatchFiles", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles", "2RatesAnalyses",
                                  true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles", "Distances",
                                  true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles", "GA", true,
                                  true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles", "GUI", true,
                                  true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles", "Miscellaneous",
                                  true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/Miscellaneous",
                                  "phylohandbook", true, true);
          Module['FS_createPath'](
              "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook",
              "datasets", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles",
                                  "ProteinAnalyses", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles", "Samplers", true,
                                  true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles",
                                  "SelectionAnalyses", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/SelectionAnalyses",
                                  "modules", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles", "TemplateModels",
                                  true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/TemplateModels",
                                  "EmpiricalAA", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/TemplateModels",
                                  "EmpiricalCodon", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles", "Utility", true,
                                  true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles", "lib", true,
                                  true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles", "libv3", true,
                                  true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/libv3",
                                  "convenience", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/libv3", "models",
                                  true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/libv3/models",
                                  "DNA", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/libv3/models",
                                  "binary", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/libv3/models",
                                  "codon", true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/libv3/models",
                                  "protein", true, true);
          Module['FS_createPath'](
              "/hyphy/TemplateBatchFiles/libv3/models/protein", "matrices",
              true, true);
          Module['FS_createPath']("/hyphy/TemplateBatchFiles/libv3", "tasks",
                                  true, true);
          Module['FS_createPath']("/hyphy", "data", true, true);
          Module['FS_createPath']("/", "tests", true, true);
          Module['FS_createPath']("/tests", "Alignment", true, true);
          Module['FS_createPath']("/tests", "Ancestors", true, true);
          Module['FS_createPath']("/tests", "BFFeatures", true, true);
          Module['FS_createPath']("/tests/BFFeatures", "Level2", true, true);
          Module['FS_createPath']("/tests", "BayesianGraphicalModels", true,
                                  true);
          Module['FS_createPath']("/tests", "HMM", true, true);
          Module['FS_createPath']("/tests", "REL", true, true);
          Module['FS_createPath']("/tests", "RegressionTesting", true, true);
          Module['FS_createPath']("/tests/RegressionTesting", "RELAX", true,
                                  true);
          Module['FS_createPath']("/tests/RegressionTesting", "res", true,
                                  true);
          Module['FS_createPath']("/tests", "Results", true, true);
          Module['FS_createPath']("/tests", "Shared", true, true);
          Module['FS_createPath']("/tests", "SimpleOptimizations", true, true);
          Module['FS_createPath']("/tests", "SpecializedOptimizations", true,
                                  true);
          Module['FS_createPath']("/tests", "Trickier", true, true);
          Module['FS_createPath']("/tests", "UnitTests", true, true);
          Module['FS_createPath']("/tests/UnitTests", "HBLCommands", true,
                                  true);
          Module['FS_createPath']("/tests/UnitTests/HBLCommands", "nested",
                                  true, true);
          Module['FS_createPath']("/tests/UnitTests/HBLCommands/nested",
                                  "nested2", true, true);
          Module['FS_createPath']("/tests/UnitTests/HBLCommands", "res", true,
                                  true);
          Module['FS_createPath']("/tests/UnitTests/HBLCommands/res", "SCFG",
                                  true, true);
          Module['FS_createPath']("/tests/UnitTests/HBLCommands", "testdata",
                                  true, true);
          Module['FS_createPath']("/tests/UnitTests/HBLCommands", "tmp", true,
                                  true);
          Module['FS_createPath']("/tests", "data", true, true);
          Module['FS_createPath']("/tests", "libv3", true, true);
          Module['FS_createPath']("/tests/libv3", "data", true, true);
          Module['FS_createPath']("/tests/libv3/data",
                                  "protgtr_fitter_alignments", true, true);
          Module['FS_createPath']("/tests/libv3", "models", true, true);
          Module['FS_createPath']("/tests/libv3", "support", true, true);

          /** @constructor */
          function DataRequest(start, end, audio) {
            this.start = start;
            this.end = end;
            this.audio = audio;
          }
          DataRequest.prototype = {
            requests : {},
            open : function(mode, name) {
              this.name = name;
              this.requests[name] = this;
              Module['addRunDependency'](`fp ${this.name}`);
            },
            send : function() {},
            onload : function() {
              var byteArray = this.byteArray.subarray(this.start, this.end);
              this.finish(byteArray);
            },
            finish : async function(byteArray) {
              var that = this;
              // canOwn this data in the filesystem, it is a slice into the heap
              // that will never change
              Module['FS_createDataFile'](this.name, null, byteArray, true,
                                          true, true);
              Module['removeRunDependency'](`fp ${that.name}`);
              this.requests[this.name] = null;
            }
          };

          var files = metadata['files'];
          for (var i = 0; i < files.length; ++i) {
            new DataRequest(files[i]['start'], files[i]['end'],
                            files[i]['audio'] || 0)
                .open('GET', files[i]['filename']);
          }

          function processPackageData(arrayBuffer) {
            assert(arrayBuffer, 'Loading data file failed.');
            assert(arrayBuffer.constructor.name === ArrayBuffer.name,
                   'bad input to processPackageData');
            var byteArray = new Uint8Array(arrayBuffer);
            var curr;
            // Reuse the bytearray from the XHR as the source for file reads.
            DataRequest.prototype.byteArray = byteArray;
            var files = metadata['files'];
            for (var i = 0; i < files.length; ++i) {
              DataRequest.prototype.requests[files[i].filename].onload();
            }
            Module['removeRunDependency']('datafile_hyphy.data');
          }
          Module['addRunDependency']('datafile_hyphy.data');

          Module['preloadResults'] ??= {};

          Module['preloadResults'][PACKAGE_NAME] = {fromCache : false};
          if (fetched) {
            processPackageData(fetched);
            fetched = null;
          } else {
            fetchedCallback = processPackageData;
          }
        }
        if (Module['calledRun']) {
          runWithFS(Module);
        } else {
          (Module['preRun'] ??= [])
              .push(runWithFS); // FS is not initialized yet, wait for it
        }
      }
      loadPackage({
        "files" : [
          {"filename" : "/hyphy/.DS_Store", "start" : 0, "end" : 10244},
          {
            "filename" : "/hyphy/GeneticCodes/Alt_Yeast_Nuclear.cod",
            "start" : 10244,
            "end" : 11052
          },
          {
            "filename" : "/hyphy/GeneticCodes/Ascidian_mtDNA.cod",
            "start" : 11052,
            "end" : 11854
          },
          {
            "filename" : "/hyphy/GeneticCodes/Blepharisma_Nuclear.cod",
            "start" : 11854,
            "end" : 12636
          },
          {
            "filename" : "/hyphy/GeneticCodes/Ciliate.cod",
            "start" : 12636,
            "end" : 13452
          },
          {
            "filename" : "/hyphy/GeneticCodes/Echinoderm_mtDNA.cod",
            "start" : 13452,
            "end" : 14256
          },
          {
            "filename" : "/hyphy/GeneticCodes/Euplotid_Nuclear.cod",
            "start" : 14256,
            "end" : 15055
          },
          {
            "filename" : "/hyphy/GeneticCodes/Flatworm_mtDNA.cod",
            "start" : 15055,
            "end" : 15855
          },
          {
            "filename" : "/hyphy/GeneticCodes/Invertebrate_mtDNA.cod",
            "start" : 15855,
            "end" : 16662
          },
          {
            "filename" : "/hyphy/GeneticCodes/Mold_mtDNA.cod",
            "start" : 16662,
            "end" : 17448
          },
          {
            "filename" : "/hyphy/GeneticCodes/Thraustochytrium_mtDNA.cod",
            "start" : 17448,
            "end" : 18266
          },
          {
            "filename" : "/hyphy/GeneticCodes/Vertebratemtdna.cod",
            "start" : 18266,
            "end" : 19609
          },
          {
            "filename" : "/hyphy/GeneticCodes/Yeast_mtDNA.cod",
            "start" : 19609,
            "end" : 20371
          },
          {"filename" : "/hyphy/README.md", "start" : 20371, "end" : 21797},
          {
            "filename" : "/hyphy/SubstitutionClasses/AAEFV/Equal",
            "start" : 21797,
            "end" : 21926
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/AAEFV/Estimated",
            "start" : 21926,
            "end" : 32336
          },
          {
            "filename" :
                "/hyphy/SubstitutionClasses/AAEFV/Observed In Data Set",
            "start" : 32336,
            "end" : 32494
          },
          {
            "filename" :
                "/hyphy/SubstitutionClasses/AAEFV/Observed In Partition",
            "start" : 32494,
            "end" : 32658
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/CodonEFV/Equal",
            "start" : 32658,
            "end" : 32960
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/CodonEFV/Observed Codon",
            "start" : 32960,
            "end" : 34794
          },
          {
            "filename" :
                "/hyphy/SubstitutionClasses/CodonEFV/Observed Nuc 3 params.",
            "start" : 34794,
            "end" : 37030
          },
          {
            "filename" :
                "/hyphy/SubstitutionClasses/CodonEFV/Observed Nuc 9 params.",
            "start" : 37030,
            "end" : 39344
          },
          {
            "filename" :
                "/hyphy/SubstitutionClasses/Heterogeneity/2 Bin Discrete",
            "start" : 39344,
            "end" : 39736
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/Heterogeneity/Beta",
            "start" : 39736,
            "end" : 40295
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/Heterogeneity/Beta-Gamma",
            "start" : 40295,
            "end" : 41278
          },
          {
            "filename" :
                "/hyphy/SubstitutionClasses/Heterogeneity/General Discrete",
            "start" : 41278,
            "end" : 43311
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/Heterogeneity/Half Normal",
            "start" : 43311,
            "end" : 44284
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/Heterogeneity/Lognormal",
            "start" : 44284,
            "end" : 44830
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/NucEFV/Equal",
            "start" : 44830,
            "end" : 44895
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/NucEFV/Estimated",
            "start" : 44895,
            "end" : 45638
          },
          {
            "filename" :
                "/hyphy/SubstitutionClasses/NucEFV/Observed In Data Set",
            "start" : 45638,
            "end" : 45795
          },
          {
            "filename" :
                "/hyphy/SubstitutionClasses/NucEFV/Observed In Partition",
            "start" : 45795,
            "end" : 45959
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/aa.bf",
            "start" : 45959,
            "end" : 48097
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/codon.bf",
            "start" : 48097,
            "end" : 54334
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/dinuc.bf",
            "start" : 54334,
            "end" : 55820
          },
          {
            "filename" : "/hyphy/SubstitutionClasses/nuc.bf",
            "start" : 55820,
            "end" : 56287
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Aminoacid/Dayhoff.mdl",
            "start" : 56287,
            "end" : 67780
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Aminoacid/EIAA.mdl",
            "start" : 67780,
            "end" : 79227
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Aminoacid/Fitness.mdl",
            "start" : 79227,
            "end" : 81413
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Aminoacid/Jones.mdl",
            "start" : 81413,
            "end" : 115749
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Aminoacid/mtREV.mdl",
            "start" : 115749,
            "end" : 181326
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Binary/F81.mdl",
            "start" : 181326,
            "end" : 182469
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Codon/GY94_3x4.mdl",
            "start" : 182469,
            "end" : 188874
          },
          {
            "filename" :
                "/hyphy/SubstitutionModels/Codon/Lineage_MG94xHKY85.mdl",
            "start" : 188874,
            "end" : 195840
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Codon/MG94.mdl",
            "start" : 195840,
            "end" : 201027
          },
          {
            "filename" :
                "/hyphy/SubstitutionModels/Codon/MG94REVOmegaCF3x4.mdl",
            "start" : 201027,
            "end" : 214920
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Codon/MG94_3x4.mdl",
            "start" : 214920,
            "end" : 220237
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Codon/MG94_HKY85x3_4.mdl",
            "start" : 220237,
            "end" : 227022
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Codon/MG94_REV_3x4.mdl",
            "start" : 227022,
            "end" : 240612
          },
          {
            "filename" :
                "/hyphy/SubstitutionModels/Codon/MG94xHKY85_3x4_2Rates.mdl",
            "start" : 240612,
            "end" : 246943
          },
          {
            "filename" :
                "/hyphy/SubstitutionModels/Codon/MG94xREV_3x4_DualRV.mdl",
            "start" : 246943,
            "end" : 258647
          },
          {
            "filename" :
                "/hyphy/SubstitutionModels/Codon/MG94xREV_3x4_DualRV_GDD.mdl",
            "start" : 258647,
            "end" : 274788
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Codon/MG94xTN93_3x4.mdl",
            "start" : 274788,
            "end" : 288298
          },
          {
            "filename" :
                "/hyphy/SubstitutionModels/Nucleotide/EFVEstimated.ibf",
            "start" : 288298,
            "end" : 289194
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Nucleotide/F81.mdl",
            "start" : 289194,
            "end" : 291102
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Nucleotide/HKY85.mdl",
            "start" : 291102,
            "end" : 293523
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Nucleotide/REV.mdl",
            "start" : 293523,
            "end" : 295993
          },
          {
            "filename" :
                "/hyphy/SubstitutionModels/Nucleotide/REVBetaGamma.mdl",
            "start" : 295993,
            "end" : 299130
          },
          {
            "filename" : "/hyphy/SubstitutionModels/Nucleotide/TrN.mdl",
            "start" : 299130,
            "end" : 301332
          },
          {
            "filename" : "/hyphy/SubstitutionModels/User/Nucleotide/JC69",
            "start" : 301332,
            "end" : 301918
          },
          {
            "filename" : "/hyphy/SubstitutionModels/User/Nucleotide/K2P",
            "start" : 301918,
            "end" : 302503
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/.DS_Store",
            "start" : 302503,
            "end" : 318891
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/2RatesAnalyses/GY94.mdl",
            "start" : 318891,
            "end" : 326941
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/MG94GY94xREV_PARRIS_syn3.mdl",
            "start" : 326941,
            "end" : 335444
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/MG94xREV.mdl",
            "start" : 335444,
            "end" : 340455
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/MG94xREVxBivariate.mdl",
            "start" : 340455,
            "end" : 344376
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/MG94xREVxBivariate_Multirate.mdl",
            "start" : 344376,
            "end" : 348984
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/PARRIS_M1.def",
            "start" : 348984,
            "end" : 349311
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/PARRIS_M2.def",
            "start" : 349311,
            "end" : 349801
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/PARRIS_M3.def",
            "start" : 349801,
            "end" : 350359
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/PARRIS_syn3.def",
            "start" : 350359,
            "end" : 351117
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/PARRIS_synvar.def",
            "start" : 351117,
            "end" : 351728
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/discreteGenerator.bf",
            "start" : 351728,
            "end" : 357047
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/discreteGeneratorNoPS.bf",
            "start" : 357047,
            "end" : 360629
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/2RatesAnalyses/gamma1.def",
            "start" : 360629,
            "end" : 361161
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/2RatesAnalyses/gamma2+Inv.def",
            "start" : 361161,
            "end" : 362090
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/2RatesAnalyses/gamma2.def",
            "start" : 362090,
            "end" : 362611
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/AAModelComparison.bf",
            "start" : 362611,
            "end" : 368075
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/AddABias.ibf",
            "start" : 368075,
            "end" : 377289
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/AnalyzeCodonData.bf",
            "start" : 377289,
            "end" : 378562
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/AnalyzeDiNucData.bf",
            "start" : 378562,
            "end" : 379223
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/AnalyzeNucDataFreq.bf",
            "start" : 379223,
            "end" : 380092
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/AnalyzeNucProtData.bf",
            "start" : 380092,
            "end" : 380960
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/BGM.bf",
            "start" : 380960,
            "end" : 401968
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/BivariateCodonRateAnalysis.bf",
            "start" : 401968,
            "end" : 405008
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/BranchSiteREL.bf",
            "start" : 405008,
            "end" : 444151
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/BranchSiteRELMultiModel.bf",
            "start" : 444151,
            "end" : 464737
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/CleanGaps.bf",
            "start" : 464737,
            "end" : 468005
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/CleanStopCodons.bf",
            "start" : 468005,
            "end" : 474404
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/ClusterAnalysis.bf",
            "start" : 474404,
            "end" : 485014
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/ClusterByDistanceRange.bf",
            "start" : 485014,
            "end" : 498805
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/CodonBivariateRateProcessor.bf",
            "start" : 498805,
            "end" : 509776
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/CodonModelCompare.bf",
            "start" : 509776,
            "end" : 535866
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/CodonToProtein.bf",
            "start" : 535866,
            "end" : 540801
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/ConvertDataFile.bf",
            "start" : 540801,
            "end" : 543086
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/DirectionalREL.bf",
            "start" : 543086,
            "end" : 561094
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/DistanceMatrix.bf",
            "start" : 561094,
            "end" : 569087
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/CodonTools.def",
            "start" : 569087,
            "end" : 569390
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/CodonTools2.def",
            "start" : 569390,
            "end" : 570029
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Distances/CodonToolsMain.def",
            "start" : 570029,
            "end" : 579844
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/JC69",
            "start" : 579844,
            "end" : 580415
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/K2P",
            "start" : 580415,
            "end" : 581424
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/K2P_RV",
            "start" : 581424,
            "end" : 582722
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Distances/Modified_Nei_Gojobori",
            "start" : 582722,
            "end" : 586475
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/Nei_Gojobori",
            "start" : 586475,
            "end" : 590027
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/PC",
            "start" : 590027,
            "end" : 590555
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/PC_MH",
            "start" : 590555,
            "end" : 591216
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/PC_RV",
            "start" : 591216,
            "end" : 592011
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/T3P",
            "start" : 592011,
            "end" : 593245
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/TN84",
            "start" : 593245,
            "end" : 594573
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/TN93",
            "start" : 594573,
            "end" : 597139
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/TN93_RV",
            "start" : 597139,
            "end" : 599062
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/Unaligned_LZ",
            "start" : 599062,
            "end" : 599841
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/Unaligned_LZ_FR",
            "start" : 599841,
            "end" : 601931
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/p_Distance",
            "start" : 601931,
            "end" : 602407
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/p_Distance_aa",
            "start" : 602407,
            "end" : 602969
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Distances/p_Distance_binary",
            "start" : 602969,
            "end" : 603393
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Distances/p_Distance_codon",
            "start" : 603393,
            "end" : 605193
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/F_ST.bf",
            "start" : 605193,
            "end" : 622949
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/FitnessAAModels.bf",
            "start" : 622949,
            "end" : 630990
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/GA/MSS-selector-codon.bf",
            "start" : 630990,
            "end" : 659648
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/GA/README.md",
            "start" : 659648,
            "end" : 664100
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/GA/processor-codon.bf",
            "start" : 664100,
            "end" : 669682
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/GARD.bf",
            "start" : 669682,
            "end" : 719954
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/GA_CHC.ibf",
            "start" : 719954,
            "end" : 731430
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/GA_CHC_Binary.ibf",
            "start" : 731430,
            "end" : 745780
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/GUI/crayolaColors.def",
            "start" : 745780,
            "end" : 745797
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Gateaux.bf",
            "start" : 745797,
            "end" : 754080
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/GateauxMR.bf",
            "start" : 754080,
            "end" : 761884
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/KHTest.bf",
            "start" : 761884,
            "end" : 768505
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/LEISR.bf",
            "start" : 768505,
            "end" : 793857
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/LHT.bf",
            "start" : 793857,
            "end" : 800376
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/LRTRecombTest.bf",
            "start" : 800376,
            "end" : 808199
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/LZ_Complexity.bf",
            "start" : 808199,
            "end" : 810190
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/LocalMolClock.bf",
            "start" : 810190,
            "end" : 815284
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/MGvsGY.bf",
            "start" : 815284,
            "end" : 831237
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/MSS-joint-fitter.bf",
            "start" : 831237,
            "end" : 846640
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/MSS-selector.bf",
            "start" : 846640,
            "end" : 875383
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/MergeSequences.bf",
            "start" : 875383,
            "end" : 876341
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/MergeSites.bf",
            "start" : 876341,
            "end" : 877307
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/BranchAPriori.bf",
            "start" : 877307,
            "end" : 884567
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/CountSubstitutions.bf",
            "start" : 884567,
            "end" : 888718
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/EffectOfTopology.bf",
            "start" : 888718,
            "end" : 896675
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/ErrorEstimates.bf",
            "start" : 896675,
            "end" : 901511
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/LRT.bf",
            "start" : 901511,
            "end" : 907567
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/LocalvsGlobal.bf",
            "start" : 907567,
            "end" : 911752
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/NeutralExpectation.bf",
            "start" : 911752,
            "end" : 915745
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/NucleotideBiases.bf",
            "start" : 915745,
            "end" : 920500
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/Phylohandbook.bf",
            "start" : 920500,
            "end" : 922795
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/WhatsInTheMean.bf",
            "start" : 922795,
            "end" : 928610
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/dSdN.bf",
            "start" : 928610,
            "end" : 937028
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/datasets/Drosophilia_adh.nex",
            "start" : 937028,
            "end" : 942014
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Miscellaneous/phylohandbook/datasets/H5N1_HA_5.nex",
            "start" : 942014,
            "end" : 951161
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/ModelTest.bf",
            "start" : 951161,
            "end" : 970728
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/MolClockAllRoots.bf",
            "start" : 970728,
            "end" : 976429
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/MolecularClock.bf",
            "start" : 976429,
            "end" : 979656
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/NeighborJoining.bf",
            "start" : 979656,
            "end" : 983816
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/NielsenYang.bf",
            "start" : 983816,
            "end" : 1016351
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/NucModelCompare.bf",
            "start" : 1016351,
            "end" : 1040976
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/PARRIS.bf",
            "start" : 1040976,
            "end" : 1083631
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/PairwiseRelativeRate.bf",
            "start" : 1083631,
            "end" : 1091073
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/PairwiseRelativeRatio.bf",
            "start" : 1091073,
            "end" : 1100919
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/PartitionDataFile.bf",
            "start" : 1100919,
            "end" : 1102713
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/ProteinAnalyses/ProteinGTRFit.bf",
            "start" : 1102713,
            "end" : 1117474
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/ProteinAnalyses/ProteinGTRFit_helper.ibf",
            "start" : 1117474,
            "end" : 1142309
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/ProteinAnalyses/README_relative_prot_rates.md",
            "start" : 1142309,
            "end" : 1144228
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/ProteinAnalyses/plusF_helper.ibf",
            "start" : 1144228,
            "end" : 1146052
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/ProteinAnalyses/relative_prot_rates.bf",
            "start" : 1146052,
            "end" : 1162003
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/ReduceDataSetMatrix.bf",
            "start" : 1162003,
            "end" : 1171110
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/RelativeRate.bf",
            "start" : 1171110,
            "end" : 1174621
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/RelativeRatio.bf",
            "start" : 1174621,
            "end" : 1180125
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SGASimProcessor.bf",
            "start" : 1180125,
            "end" : 1193375
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SGEmulator.bf",
            "start" : 1193375,
            "end" : 1216777
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SGEmulator_MF.bf",
            "start" : 1216777,
            "end" : 1236307
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SGIvL.bf",
            "start" : 1236307,
            "end" : 1255163
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SM-2019.bf",
            "start" : 1255163,
            "end" : 1267842
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Samplers/lhc-ErrorEst.bf",
            "start" : 1267842,
            "end" : 1269440
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Samplers/lhc.bf",
            "start" : 1269440,
            "end" : 1275963
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Samplers/lhc_supp.ibf",
            "start" : 1275963,
            "end" : 1277200
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Samplers/sir.bf",
            "start" : 1277200,
            "end" : 1280467
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Samplers/srs-ErrorEst.ibf",
            "start" : 1280467,
            "end" : 1282847
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Samplers/srs.ibf",
            "start" : 1282847,
            "end" : 1291863
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SandNSAmbigs.bf",
            "start" : 1291863,
            "end" : 1301390
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/.DS_Store",
            "start" : 1301390,
            "end" : 1309586
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/BUSTED-PH.bf",
            "start" : 1309586,
            "end" : 1329193
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf",
            "start" : 1329193,
            "end" : 1420401
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/BranchSiteREL.bf",
            "start" : 1420401,
            "end" : 1463239
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SelectionAnalyses/FADE.bf",
            "start" : 1463239,
            "end" : 1509914
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SelectionAnalyses/FEL.bf",
            "start" : 1509914,
            "end" : 1572259
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SelectionAnalyses/FUBAR.bf",
            "start" : 1572259,
            "end" : 1613535
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/FitMultiModel.bf",
            "start" : 1613535,
            "end" : 1639908
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SelectionAnalyses/MEME.bf",
            "start" : 1639908,
            "end" : 1721555
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SelectionAnalyses/PRIME.bf",
            "start" : 1721555,
            "end" : 1775107
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/RELAX-Groups.bf",
            "start" : 1775107,
            "end" : 1793055
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SelectionAnalyses/RELAX.bf",
            "start" : 1793055,
            "end" : 1892760
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SelectionAnalyses/SLAC.bf",
            "start" : 1892760,
            "end" : 1927242
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/SingleOmega.bf",
            "start" : 1927242,
            "end" : 1934093
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/aBSREL.bf",
            "start" : 1934093,
            "end" : 2003608
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/contrast-fel.bf",
            "start" : 2003608,
            "end" : 2049575
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/contrast-meme.bf",
            "start" : 2049575,
            "end" : 2099776
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/error-filter.bf",
            "start" : 2099776,
            "end" : 2112872
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/grid.json",
            "start" : 2112872,
            "end" : 2342212
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/modules/grid_compute.ibf",
            "start" : 2342212,
            "end" : 2361873
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/modules/io_functions.ibf",
            "start" : 2361873,
            "end" : 2382110
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/modules/selection_lib.ibf",
            "start" : 2382110,
            "end" : 2386870
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SelectionAnalyses/modules/shared-load-file.bf",
            "start" : 2386870,
            "end" : 2424362
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SelectionLRT.bf",
            "start" : 2424362,
            "end" : 2435213
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SeqAlignShared.ibf",
            "start" : 2435213,
            "end" : 2520237
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SeqAlignmentCodon.bf",
            "start" : 2520237,
            "end" : 2520406
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/SeqAlignmentCodonShared.ibf",
            "start" : 2520406,
            "end" : 2528302
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SeqAlignmentNucShared.ibf",
            "start" : 2528302,
            "end" : 2535515
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SequentialAddition.bf",
            "start" : 2535515,
            "end" : 2541054
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SequentialAddition.ibf",
            "start" : 2541054,
            "end" : 2556322
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SimilarityPlot.bf",
            "start" : 2556322,
            "end" : 2560956
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SimmondsAI.bf",
            "start" : 2560956,
            "end" : 2568755
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SingleBreakpointRecomb.bf",
            "start" : 2568755,
            "end" : 2598188
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/SlidingWindowAnalysis.bf",
            "start" : 2598188,
            "end" : 2604733
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/.DS_Store",
            "start" : 2604733,
            "end" : 2612929
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/190",
            "start" : 2612929,
            "end" : 2614334
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/BranchSiteTemplate.mdl",
            "start" : 2614334,
            "end" : 2621453
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/CF3x4.bf",
            "start" : 2621453,
            "end" : 2624188
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/Custom_AA.mdl",
            "start" : 2624188,
            "end" : 2627710
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/Custom_AA_empirical.mdl",
            "start" : 2627710,
            "end" : 2632067
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/Dayhoff.mdl",
            "start" : 2632067,
            "end" : 2666595
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/Dayhoff_F.mdl",
            "start" : 2666595,
            "end" : 2700638
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/Default",
            "start" : 2700638,
            "end" : 2701505
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/ECM+F+omega.mdl",
            "start" : 2701505,
            "end" : 2705660
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/ECM+F.mdl",
            "start" : 2705660,
            "end" : 2709813
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/ECM+F3x4.mdl",
            "start" : 2709813,
            "end" : 2714447
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/ECM+MLFREQS.mdl",
            "start" : 2714447,
            "end" : 2718752
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/ECM.mdl",
            "start" : 2718752,
            "end" : 2722906
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/EIAA.mdl",
            "start" : 2722906,
            "end" : 2724542
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EIAAFreq.mdl",
            "start" : 2724542,
            "end" : 2726782
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/EX.dat",
            "start" : 2726782,
            "end" : 2729509
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/EX.mdl",
            "start" : 2729509,
            "end" : 2764181
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/BLOSUM62",
            "start" : 2764181,
            "end" : 2770589
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/Dayhoff",
            "start" : 2770589,
            "end" : 2776999
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/H5N1",
            "start" : 2776999,
            "end" : 2785009
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/HIVBetween",
            "start" : 2785009,
            "end" : 2788516
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/HIVWithin",
            "start" : 2788516,
            "end" : 2791749
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/IAV",
            "start" : 2791749,
            "end" : 2798159
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/JTT",
            "start" : 2798159,
            "end" : 2806250
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/LCAP",
            "start" : 2806250,
            "end" : 2842890
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/LG",
            "start" : 2842890,
            "end" : 2849298
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/MtArt",
            "start" : 2849298,
            "end" : 2855706
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/VT",
            "start" : 2855706,
            "end" : 2862114
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/WAG",
            "start" : 2862114,
            "end" : 2868524
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/cpREV",
            "start" : 2868524,
            "end" : 2874932
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/modellist.ibf",
            "start" : 2874932,
            "end" : 2885079
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/mtMAM",
            "start" : 2885079,
            "end" : 2891489
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/mtREV24",
            "start" : 2891489,
            "end" : 2897899
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalAA/rtREV",
            "start" : 2897899,
            "end" : 2904309
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalCodon/KHG_ECM",
            "start" : 2904309,
            "end" : 2942286
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/EmpiricalCodon/KHG_ECMu",
            "start" : 2942286,
            "end" : 2980315
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/F81.mdl",
            "start" : 2980315,
            "end" : 2981574
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/F81_binary.mdl",
            "start" : 2981574,
            "end" : 2982765
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/F84.mdl",
            "start" : 2982765,
            "end" : 2985151
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/F84P.mdl",
            "start" : 2985151,
            "end" : 2988906
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/GRM.mdl",
            "start" : 2988906,
            "end" : 2990578
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/GY94.ibf",
            "start" : 2990578,
            "end" : 2995563
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/GY94.mdl",
            "start" : 2995563,
            "end" : 2997436
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/GY94customF3x4.mdl",
            "start" : 2997436,
            "end" : 3004613
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/GY94w9.mdl",
            "start" : 3004613,
            "end" : 3006792
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/HIVbetween+F.mdl",
            "start" : 3006792,
            "end" : 3007147
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/HIVbetween.ibf",
            "start" : 3007147,
            "end" : 3011096
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/HIVbetween.mdl",
            "start" : 3011096,
            "end" : 3011693
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/HIVwithin+F.mdl",
            "start" : 3011693,
            "end" : 3012044
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/HIVwithin.ibf",
            "start" : 3012044,
            "end" : 3016116
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/HIVwithin.mdl",
            "start" : 3016116,
            "end" : 3016606
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/HKY85.mdl",
            "start" : 3016606,
            "end" : 3018135
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/JC69.mdl",
            "start" : 3018135,
            "end" : 3019391
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/JC69_binary.mdl",
            "start" : 3019391,
            "end" : 3020566
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/Jones.mdl",
            "start" : 3020566,
            "end" : 3055653
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/Jones_F.mdl",
            "start" : 3055653,
            "end" : 3089662
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/K2P.mdl",
            "start" : 3089662,
            "end" : 3091084
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/K3ST.mdl",
            "start" : 3091084,
            "end" : 3092597
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/LCAP.mdl",
            "start" : 3092597,
            "end" : 3101768
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/MEC.mdl",
            "start" : 3101768,
            "end" : 3109781
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/MG94.mdl",
            "start" : 3109781,
            "end" : 3115635
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/MG94custom.mdl",
            "start" : 3115635,
            "end" : 3125207
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/MG94customCF3x4.mdl",
            "start" : 3125207,
            "end" : 3134978
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/MG94customF1x4.mdl",
            "start" : 3134978,
            "end" : 3144560
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/MG94customFreqs.mdl",
            "start" : 3144560,
            "end" : 3153627
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/MG94wAA.mdl",
            "start" : 3153627,
            "end" : 3159709
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/MG94wAAF61.mdl",
            "start" : 3159709,
            "end" : 3165062
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/MG94wAAF61multiple.mdl",
            "start" : 3165062,
            "end" : 3170646
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/MG94wAAFreqs.mdl",
            "start" : 3170646,
            "end" : 3176145
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/MG94wAAUserFreqs.mdl",
            "start" : 3176145,
            "end" : 3181554
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/MG94wEX.mdl",
            "start" : 3181554,
            "end" : 3187050
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/MG94with9freqs.mdl",
            "start" : 3187050,
            "end" : 3193491
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/MG94x2.mdl",
            "start" : 3193491,
            "end" : 3204553
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/MGFreqsEstimator.ibf",
            "start" : 3204553,
            "end" : 3206086
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/MGwAA.ibf",
            "start" : 3206086,
            "end" : 3209442
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/MGwEX.ibf",
            "start" : 3209442,
            "end" : 3212130
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/NRM+Freqs.mdl",
            "start" : 3212130,
            "end" : 3213714
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/NRM.mdl",
            "start" : 3213714,
            "end" : 3215453
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/PCH",
            "start" : 3215453,
            "end" : 3216320
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/RNA16.mdl",
            "start" : 3216320,
            "end" : 3222148
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/RNA16A",
            "start" : 3222148,
            "end" : 3223628
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/RNA16SvH.mdl",
            "start" : 3223628,
            "end" : 3227115
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/RNAEqualInput",
            "start" : 3227115,
            "end" : 3227664
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/RNAF81",
            "start" : 3227664,
            "end" : 3228580
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/RNAMuse95.mdl",
            "start" : 3228580,
            "end" : 3232329
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/RNAREV",
            "start" : 3232329,
            "end" : 3233138
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/RNAREV_1",
            "start" : 3233138,
            "end" : 3233923
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/S4",
            "start" : 3233923,
            "end" : 3234789
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/STGRM.mdl",
            "start" : 3234789,
            "end" : 3236445
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/TrN.mdl",
            "start" : 3236445,
            "end" : 3238009
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/UniversalCode.def",
            "start" : 3238009,
            "end" : 3243128
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/WAG.mdl",
            "start" : 3243128,
            "end" : 3271696
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/WAG_F.mdl",
            "start" : 3271696,
            "end" : 3299819
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/Yang2000Distributions.def",
            "start" : 3299819,
            "end" : 3310496
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/chooseGeneticCode.def",
            "start" : 3310496,
            "end" : 3342638
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/custm4x4.mdl",
            "start" : 3342638,
            "end" : 3352009
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/defineGamma.mdl",
            "start" : 3352009,
            "end" : 3362990
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/defineHM.mdl",
            "start" : 3362990,
            "end" : 3373187
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/fitness.mdl",
            "start" : 3373187,
            "end" : 3375065
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/modelParameters.mdl",
            "start" : 3375065,
            "end" : 3375894
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/modelParameters2.mdl",
            "start" : 3375894,
            "end" : 3376440
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/modelParameters3.mdl",
            "start" : 3376440,
            "end" : 3377099
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/modelParameters4.mdl",
            "start" : 3377099,
            "end" : 3378798
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/modelParameters5.mdl",
            "start" : 3378798,
            "end" : 3379458
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/models.lst",
            "start" : 3379458,
            "end" : 3390444
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/mtMAM.mdl",
            "start" : 3390444,
            "end" : 3406822
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/mtMAM_F.mdl",
            "start" : 3406822,
            "end" : 3422819
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/mtREV.mdl",
            "start" : 3422819,
            "end" : 3424293
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/mtREV_24.mdl",
            "start" : 3424293,
            "end" : 3452517
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/mtREV_24_F.mdl",
            "start" : 3452517,
            "end" : 3480374
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TemplateModels/reducedREV.mdl",
            "start" : 3480374,
            "end" : 3483434
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/rtREV.mdl",
            "start" : 3483434,
            "end" : 3512168
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TemplateModels/rtREV_F.mdl",
            "start" : 3512168,
            "end" : 3540514
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TestBranchDNDS.bf",
            "start" : 3540514,
            "end" : 3552788
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TestClade.bf",
            "start" : 3552788,
            "end" : 3576502
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TestCladeMeans.bf",
            "start" : 3576502,
            "end" : 3586797
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TopologySearch.bf",
            "start" : 3586797,
            "end" : 3597580
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TopologySearchConstrained.bf",
            "start" : 3597580,
            "end" : 3606670
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/TreeCorrelationCoefficients.bf",
            "start" : 3606670,
            "end" : 3614639
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/TreeTools.ibf",
            "start" : 3614639,
            "end" : 3635317
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/UpperBound.bf",
            "start" : 3635317,
            "end" : 3636962
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Utility/AncestralMapper.bf",
            "start" : 3636962,
            "end" : 3660577
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Utility/BranchLengthFitters.bf",
            "start" : 3660577,
            "end" : 3663293
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Utility/CoalescentPostProcessor.bf",
            "start" : 3663293,
            "end" : 3675664
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Utility/CodonTools.bf",
            "start" : 3675664,
            "end" : 3677331
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Utility/DescriptiveStatistics.bf",
            "start" : 3677331,
            "end" : 3680104
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Utility/GrabBag.bf",
            "start" : 3680104,
            "end" : 3702641
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Utility/HXB2Mapper.bf",
            "start" : 3702641,
            "end" : 3723018
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Utility/LocalMGREV.bf",
            "start" : 3723018,
            "end" : 3726208
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Utility/LocalMGREVMLFreqs.bf",
            "start" : 3726208,
            "end" : 3731557
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Utility/MPITools.bf",
            "start" : 3731557,
            "end" : 3734174
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Utility/NJ.bf",
            "start" : 3734174,
            "end" : 3740259
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Utility/PS_Plotters.bf",
            "start" : 3740259,
            "end" : 3791312
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Utility/PostScript.bf",
            "start" : 3791312,
            "end" : 3803511
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Utility/ProbabilityDistributions.bf",
            "start" : 3803511,
            "end" : 3805325
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Utility/ReadDelimitedFiles.bf",
            "start" : 3805325,
            "end" : 3816910
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/Utility/TreeFunctions.bf",
            "start" : 3816910,
            "end" : 3818086
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/Utility/WriteDelimitedFiles.bf",
            "start" : 3818086,
            "end" : 3819280
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/WANC.bf",
            "start" : 3819280,
            "end" : 3828596
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/YangNielsenBranchSite2005.bf",
            "start" : 3828596,
            "end" : 3843251
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/_BMS_Aux.ibf",
            "start" : 3843251,
            "end" : 3847780
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/_CMS_Aux.ibf",
            "start" : 3847780,
            "end" : 3856635
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/_MFReader_.ibf",
            "start" : 3856635,
            "end" : 3863634
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/_tipDater.ibf",
            "start" : 3863634,
            "end" : 3865528
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/bayesgraph.ibf",
            "start" : 3865528,
            "end" : 3885972
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/binomial.ibf",
            "start" : 3885972,
            "end" : 3887619
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/branchSwappingFunctions.bf",
            "start" : 3887619,
            "end" : 3896573
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/categoryEcho.bf",
            "start" : 3896573,
            "end" : 3901014
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/chooseDistanceFormula.def",
            "start" : 3901014,
            "end" : 3905039
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/dNdSBivariateRateAnalysis.bf",
            "start" : 3905039,
            "end" : 3927059
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/dNdSDistributionComparison.bf",
            "start" : 3927059,
            "end" : 3956371
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/dNdSRateAnalysis.bf",
            "start" : 3956371,
            "end" : 3983162
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/dNdSResultProcessor.bf",
            "start" : 3983162,
            "end" : 4016031
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/dSdNTreeTools.ibf",
            "start" : 4016031,
            "end" : 4018692
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/distanceMethodNPBootstrap.bf",
            "start" : 4018692,
            "end" : 4028388
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/distanceRMethodNPBootstrap.bf",
            "start" : 4028388,
            "end" : 4036247
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/doNNISwap.bf",
            "start" : 4036247,
            "end" : 4039328
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/doSPRSwap.bf",
            "start" : 4039328,
            "end" : 4042173
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/files.lst",
            "start" : 4042173,
            "end" : 4055572
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/globalChecker.ibf",
            "start" : 4055572,
            "end" : 4056083
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/heuristicMethodNPBootstrap.bf",
            "start" : 4056083,
            "end" : 4061976
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/last.date",
            "start" : 4061976,
            "end" : 4061991
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/lib/label-tree.bf",
            "start" : 4061991,
            "end" : 4069131
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/lib/remove-duplicates.bf",
            "start" : 4069131,
            "end" : 4073974
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/lib/trim-label-tree.bf",
            "start" : 4073974,
            "end" : 4083592
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/lib/trim-tree.bf",
            "start" : 4083592,
            "end" : 4090612
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/.DS_Store",
            "start" : 4090612,
            "end" : 4104952
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/IOFunctions.bf",
            "start" : 4104952,
            "end" : 4129523
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/UtilityFunctions.bf",
            "start" : 4129523,
            "end" : 4165316
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/all-terms.bf",
            "start" : 4165316,
            "end" : 4185915
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/convenience/math.bf",
            "start" : 4185915,
            "end" : 4193621
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/convenience/matrix.bf",
            "start" : 4193621,
            "end" : 4194100
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/convenience/random.bf",
            "start" : 4194100,
            "end" : 4206570
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/convenience/regexp.bf",
            "start" : 4206570,
            "end" : 4212261
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/function-loader.bf",
            "start" : 4212261,
            "end" : 4212972
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/models/.DS_Store",
            "start" : 4212972,
            "end" : 4223216
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/models/DNA.bf",
            "start" : 4223216,
            "end" : 4226939
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/DNA/GTR-MutSel.bf",
            "start" : 4226939,
            "end" : 4230330
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/models/DNA/GTR.bf",
            "start" : 4230330,
            "end" : 4233787
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/models/DNA/HKY85.bf",
            "start" : 4233787,
            "end" : 4237638
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/models/DNA/JC69.bf",
            "start" : 4237638,
            "end" : 4240604
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/models/binary.bf",
            "start" : 4240604,
            "end" : 4243986
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/binary/charbinary.bf",
            "start" : 4243986,
            "end" : 4247504
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/binary/empirical.bf",
            "start" : 4247504,
            "end" : 4248611
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/models/codon.bf",
            "start" : 4248611,
            "end" : 4256755
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/codon/.DS_Store",
            "start" : 4256755,
            "end" : 4262903
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/codon/BS_REL.bf",
            "start" : 4262903,
            "end" : 4293229
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/codon/MG_GTR.bf",
            "start" : 4293229,
            "end" : 4298023
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/codon/MG_REV.bf",
            "start" : 4298023,
            "end" : 4313357
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/codon/MG_REV_MH.bf",
            "start" : 4313357,
            "end" : 4317936
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/codon/MG_REV_PROPERTIES.bf",
            "start" : 4317936,
            "end" : 4347196
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/codon/MG_REV_PROPERTIES_BSREL.bf",
            "start" : 4347196,
            "end" : 4353222
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/codon/MG_REV_TRIP.bf",
            "start" : 4353222,
            "end" : 4359143
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/models/codon/MSS.bf",
            "start" : 4359143,
            "end" : 4382691
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/frequencies.bf",
            "start" : 4382691,
            "end" : 4407191
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/model_functions.bf",
            "start" : 4407191,
            "end" : 4438132
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/models/parameters.bf",
            "start" : 4438132,
            "end" : 4470683
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/models/protein.bf",
            "start" : 4470683,
            "end" : 4476178
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/.DS_Store",
            "start" : 4476178,
            "end" : 4482326
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/REV.bf",
            "start" : 4482326,
            "end" : 4489786
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/empirical.bf",
            "start" : 4489786,
            "end" : 4521885
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/matrices/HIV.ibf",
            "start" : 4521885,
            "end" : 4535127
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/matrices/JC69.ibf",
            "start" : 4535127,
            "end" : 4541914
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/matrices/JTT.ibf",
            "start" : 4541914,
            "end" : 4548910
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/matrices/LG.ibf",
            "start" : 4548910,
            "end" : 4556662
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/matrices/WAG.ibf",
            "start" : 4556662,
            "end" : 4564584
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/matrices/gcpREV.ibf",
            "start" : 4564584,
            "end" : 4571804
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/matrices/mt.ibf",
            "start" : 4571804,
            "end" : 4597133
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/protein/matrix2dict_aux.bf",
            "start" : 4597133,
            "end" : 4623023
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/models/rate_variation.bf",
            "start" : 4623023,
            "end" : 4635524
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/stats.bf",
            "start" : 4635524,
            "end" : 4636534
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/tasks/.DS_Store",
            "start" : 4636534,
            "end" : 4642682
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/tasks/alignments.bf",
            "start" : 4642682,
            "end" : 4688147
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/tasks/ancestral.bf",
            "start" : 4688147,
            "end" : 4735731
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/tasks/bayesgraph.ibf",
            "start" : 4735731,
            "end" : 4763461
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/tasks/distances.bf",
            "start" : 4763461,
            "end" : 4767346
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/tasks/estimators.bf",
            "start" : 4767346,
            "end" : 4823613
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/libv3/tasks/genetic_code.bf",
            "start" : 4823613,
            "end" : 4839745
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/tasks/mapping.bf",
            "start" : 4839745,
            "end" : 4863157
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/tasks/mpi.bf",
            "start" : 4863157,
            "end" : 4880778
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/libv3/tasks/trees.bf",
            "start" : 4880778,
            "end" : 4928438
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/molclockBootstrap.bf",
            "start" : 4928438,
            "end" : 4934433
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/molerate.bf",
            "start" : 4934433,
            "end" : 4968837
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/pairwiseDistanceEstimator.ibf",
            "start" : 4968837,
            "end" : 4979361
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/pairwiseDistanceEstimatorCounter.ibf",
            "start" : 4979361,
            "end" : 4979541
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/partitionSequences.ibf",
            "start" : 4979541,
            "end" : 4981198
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/posteriors.ibf",
            "start" : 4981198,
            "end" : 4981894
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/qndhelper1.ibf",
            "start" : 4981894,
            "end" : 5017989
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/qndhelper1_mf.ibf",
            "start" : 5017989,
            "end" : 5025612
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/qndhelper2.ibf",
            "start" : 5025612,
            "end" : 5039090
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/qndhelper2_mf.ibf",
            "start" : 5039090,
            "end" : 5053236
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/qndhelper3.ibf",
            "start" : 5053236,
            "end" : 5057645
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/qndhelper4.ibf",
            "start" : 5057645,
            "end" : 5061998
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/queryTree.bf",
            "start" : 5061998,
            "end" : 5064832
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/readIndexFile.bf",
            "start" : 5064832,
            "end" : 5065280
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/relative_nucleotide_rates.bf",
            "start" : 5065280,
            "end" : 5077908
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/relrateBootstrap.bf",
            "start" : 5077908,
            "end" : 5088080
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/relratioBootstrap.bf",
            "start" : 5088080,
            "end" : 5092582
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/selectModelParameters.bf",
            "start" : 5092582,
            "end" : 5093031
          },
          {
            "filename" : "/hyphy/TemplateBatchFiles/simpleBootstrap.bf",
            "start" : 5093031,
            "end" : 5102988
          },
          {
            "filename" :
                "/hyphy/TemplateBatchFiles/temporal-summary-sites-clades.py",
            "start" : 5102988,
            "end" : 5105113
          },
          {
            "filename" : "/hyphy/data/brown.nuc",
            "start" : 5105113,
            "end" : 5109793
          },
          {
            "filename" : "/hyphy/data/integrase_BDA.nex",
            "start" : 5109793,
            "end" : 5120398
          },
          {
            "filename" : "/hyphy/data/p51.aa",
            "start" : 5120398,
            "end" : 5124620
          },
          {
            "filename" : "/hyphy/data/p51.nex",
            "start" : 5124620,
            "end" : 5135796
          },
          {"filename" : "/tests/.DS_Store", "start" : 5135796, "end" : 5156280},
          {
            "filename" : "/tests/Alignment/codon-crash.bf",
            "start" : 5156280,
            "end" : 7903730
          },
          {
            "filename" : "/tests/Alignment/codon-test.bf",
            "start" : 7903730,
            "end" : 7907915
          },
          {
            "filename" : "/tests/Alignment/nuc-test.bf",
            "start" : 7907915,
            "end" : 7912346
          },
          {
            "filename" : "/tests/Ancestors/.DS_Store",
            "start" : 7912346,
            "end" : 7918494
          },
          {
            "filename" : "/tests/Ancestors/CodonAncestors.bf",
            "start" : 7918494,
            "end" : 7954397
          },
          {
            "filename" : "/tests/Ancestors/LeafProbs.bf",
            "start" : 7954397,
            "end" : 7968581
          },
          {
            "filename" : "/tests/Ancestors/NucAncestors.bf",
            "start" : 7968581,
            "end" : 7993790
          },
          {
            "filename" : "/tests/Ancestors/NucAncestorsMultiplePartitions.bf",
            "start" : 7993790,
            "end" : 8082078
          },
          {
            "filename" : "/tests/Ancestors/NucRVAncestors.bf",
            "start" : 8082078,
            "end" : 8097004
          },
          {
            "filename" : "/tests/BFFeatures/Level1.bf",
            "start" : 8097004,
            "end" : 8097037
          },
          {
            "filename" : "/tests/BFFeatures/Level2/File1.bf",
            "start" : 8097037,
            "end" : 8097063
          },
          {
            "filename" : "/tests/BFFeatures/Level2/File2.bf",
            "start" : 8097063,
            "end" : 8097109
          },
          {
            "filename" : "/tests/BFFeatures/TreeSplits.bf",
            "start" : 8097109,
            "end" : 8098517
          },
          {
            "filename" : "/tests/BFFeatures/kernel.dump",
            "start" : 8098517,
            "end" : 9820901
          },
          {
            "filename" : "/tests/BFFeatures/simplify.bf",
            "start" : 9820901,
            "end" : 9821749
          },
          {
            "filename" : "/tests/BayesianGraphicalModels/TestBGM.bf",
            "start" : 9821749,
            "end" : 9825925
          },
          {
            "filename" : "/tests/BayesianGraphicalModels/alarm.xml",
            "start" : 9825925,
            "end" : 9841281
          },
          {
            "filename" : "/tests/HMM/RateHMM.bf",
            "start" : 9841281,
            "end" : 9845090
          },
          {
            "filename" : "/tests/HMM/SmallNuc.bf",
            "start" : 9845090,
            "end" : 9859609
          },
          {
            "filename" : "/tests/HMM/TreeHMM.bf",
            "start" : 9859609,
            "end" : 9863042
          },
          {
            "filename" : "/tests/REL/BS-REL.bf",
            "start" : 9863042,
            "end" : 9864214
          },
          {
            "filename" : "/tests/REL/BUSTED.bf",
            "start" : 9864214,
            "end" : 9865278
          },
          {
            "filename" : "/tests/REL/GTR_G_I.bf",
            "start" : 9865278,
            "end" : 9883832
          },
          {
            "filename" : "/tests/REL/IntermediateNucRel.bf",
            "start" : 9883832,
            "end" : 9958791
          },
          {
            "filename" : "/tests/REL/ModelMixture.bf",
            "start" : 9958791,
            "end" : 10135261
          },
          {
            "filename" : "/tests/REL/MultiplePartitions.bf",
            "start" : 10135261,
            "end" : 10438214
          },
          {
            "filename" : "/tests/REL/NY.bf",
            "start" : 10438214,
            "end" : 10439735
          },
          {
            "filename" : "/tests/REL/SmallNucRel.bf",
            "start" : 10439735,
            "end" : 10582845
          },
          {
            "filename" : "/tests/RegressionTesting/.DS_Store",
            "start" : 10582845,
            "end" : 10588993
          },
          {
            "filename" : "/tests/RegressionTesting/60901d3d9e9f70521111dc5e",
            "start" : 10588993,
            "end" : 10606939
          },
          {
            "filename" : "/tests/RegressionTesting/ClearConstraintsBug.bf",
            "start" : 10606939,
            "end" : 10658785
          },
          {
            "filename" : "/tests/RegressionTesting/LFReuse.bf",
            "start" : 10658785,
            "end" : 10694891
          },
          {
            "filename" : "/tests/RegressionTesting/LocalReferenceBug.bf",
            "start" : 10694891,
            "end" : 10695089
          },
          {
            "filename" : "/tests/RegressionTesting/MatrixElementExpBug.bf",
            "start" : 10695089,
            "end" : 10695251
          },
          {
            "filename" : "/tests/RegressionTesting/ParseNexus.bf",
            "start" : 10695251,
            "end" : 10696798
          },
          {
            "filename" : "/tests/RegressionTesting/RELAX/segfault.nex",
            "start" : 10696798,
            "end" : 10739105
          },
          {
            "filename" : "/tests/RegressionTesting/RELAX/segfault.nex.tre",
            "start" : 10739105,
            "end" : 10742841
          },
          {
            "filename" : "/tests/RegressionTesting/RELAX/wrapper.bf",
            "start" : 10742841,
            "end" : 10743125
          },
          {
            "filename" : "/tests/RegressionTesting/expModelCrash.bf",
            "start" : 10743125,
            "end" : 10820817
          },
          {
            "filename" : "/tests/RegressionTesting/formula-to-str-crash.bf",
            "start" : 10820817,
            "end" : 10821037
          },
          {
            "filename" : "/tests/RegressionTesting/global_export_miss.bf",
            "start" : 10821037,
            "end" : 10984790
          },
          {
            "filename" : "/tests/RegressionTesting/nan-exp.bf",
            "start" : 10984790,
            "end" : 11117662
          },
          {
            "filename" : "/tests/RegressionTesting/res/69genes.test.nex",
            "start" : 11117662,
            "end" : 11231037
          },
          {
            "filename" : "/tests/RegressionTesting/unaryMinus.bf",
            "start" : 11231037,
            "end" : 11231145
          },
          {
            "filename" : "/tests/RegressionTesting/update_listener.bf",
            "start" : 11231145,
            "end" : 11334843
          },
          {
            "filename" : "/tests/Results/HIVSweden.out",
            "start" : 11334843,
            "end" : 11348158
          },
          {
            "filename" : "/tests/Results/HIVSweden.out_MODEL_-1.nex",
            "start" : 11348158,
            "end" : 11370934
          },
          {
            "filename" : "/tests/Results/HIVSweden.out_MODEL_0.nex",
            "start" : 11370934,
            "end" : 11393911
          },
          {
            "filename" : "/tests/Results/HIVSweden.out_MODEL_1.nex",
            "start" : 11393911,
            "end" : 11417034
          },
          {
            "filename" : "/tests/Shared/REL_utils.bf",
            "start" : 11417034,
            "end" : 11417530
          },
          {
            "filename" : "/tests/Shared/TestInstrumentation.bf",
            "start" : 11417530,
            "end" : 11420868
          },
          {
            "filename" : "/tests/SimpleOptimizations/.DS_Store",
            "start" : 11420868,
            "end" : 11427016
          },
          {
            "filename" : "/tests/SimpleOptimizations/Deuterostomes.nex",
            "start" : 11427016,
            "end" : 11455526
          },
          {
            "filename" : "/tests/SimpleOptimizations/IntermediateCodon.bf",
            "start" : 11455526,
            "end" : 11506222
          },
          {
            "filename" : "/tests/SimpleOptimizations/IntermediateNuc.bf",
            "start" : 11506222,
            "end" : 11509502
          },
          {
            "filename" : "/tests/SimpleOptimizations/IntermediateProtein.bf",
            "start" : 11509502,
            "end" : 11692305
          },
          {
            "filename" : "/tests/SimpleOptimizations/LargeNuc.bf",
            "start" : 11692305,
            "end" : 14990828
          },
          {
            "filename" : "/tests/SimpleOptimizations/SmallCodon.bf",
            "start" : 14990828,
            "end" : 15026882
          },
          {
            "filename" : "/tests/SimpleOptimizations/SmallCodonLocal.bf",
            "start" : 15026882,
            "end" : 15063204
          },
          {
            "filename" : "/tests/SimpleOptimizations/TwoSequenceTest.bf",
            "start" : 15063204,
            "end" : 15065415
          },
          {
            "filename" :
                "/tests/SpecializedOptimizations/AminoAcidPartitions.bf",
            "start" : 15065415,
            "end" : 15214962
          },
          {
            "filename" : "/tests/SpecializedOptimizations/MEME.bf",
            "start" : 15214962,
            "end" : 15536060
          },
          {
            "filename" :
                "/tests/SpecializedOptimizations/SingleSiteTemplate.bf",
            "start" : 15536060,
            "end" : 15890892
          },
          {
            "filename" : "/tests/SpecializedOptimizations/SiteLikelihood.bf",
            "start" : 15890892,
            "end" : 15958880
          },
          {
            "filename" :
                "/tests/Trickier/NonReversibleWithDynamicFrequencies.bf",
            "start" : 15958880,
            "end" : 15959778
          },
          {
            "filename" : "/tests/UnitTests/.DS_Store",
            "start" : 15959778,
            "end" : 15967974
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Abs.bf",
            "start" : 15967974,
            "end" : 15969118
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Add.bf",
            "start" : 15969118,
            "end" : 15975555
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Ampersand.bf",
            "start" : 15975555,
            "end" : 15976883
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Arctan.bf",
            "start" : 15976883,
            "end" : 15979248
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Asterisk.bf",
            "start" : 15979248,
            "end" : 15982680
          },
          {
            "filename" :
                "/tests/UnitTests/HBLCommands/BayesianGraphicalModel.bf",
            "start" : 15982680,
            "end" : 15985341
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Beta.bf",
            "start" : 15985341,
            "end" : 15987801
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/BranchLength.bf",
            "start" : 15987801,
            "end" : 15993444
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Branchcount.bf",
            "start" : 15993444,
            "end" : 15995897
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Branchname.bf",
            "start" : 15995897,
            "end" : 15999518
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/CChi2.bf",
            "start" : 15999518,
            "end" : 16003348
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/CGammaDist.bf",
            "start" : 16003348,
            "end" : 16005706
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Call.bf",
            "start" : 16005706,
            "end" : 16008433
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Caret.bf",
            "start" : 16008433,
            "end" : 16013072
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Category.bf",
            "start" : 16013072,
            "end" : 16013848
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Columns.bf",
            "start" : 16013848,
            "end" : 16016226
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Continue.bf",
            "start" : 16016226,
            "end" : 16018327
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Cos.bf",
            "start" : 16018327,
            "end" : 16020629
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/DataSet.bf",
            "start" : 16020629,
            "end" : 16023833
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/DataSetFilter.bf",
            "start" : 16023833,
            "end" : 16026637
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/DeleteObject.bf",
            "start" : 16026637,
            "end" : 16027325
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Differentiate.bf",
            "start" : 16027325,
            "end" : 16036877
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Divide.bf",
            "start" : 16036877,
            "end" : 16040343
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Do.bf",
            "start" : 16040343,
            "end" : 16041608
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Dollarsign.bf",
            "start" : 16041608,
            "end" : 16044260
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/DoubleAmpersand.bf",
            "start" : 16044260,
            "end" : 16047636
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/DoubleEquals.bf",
            "start" : 16047636,
            "end" : 16050270
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/DoubleVerticalBar.bf",
            "start" : 16050270,
            "end" : 16052153
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Eigensystem.bf",
            "start" : 16052153,
            "end" : 16054595
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Erf.bf",
            "start" : 16054595,
            "end" : 16057329
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Eval.bf",
            "start" : 16057329,
            "end" : 16059259
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/ExclamationPoint.bf",
            "start" : 16059259,
            "end" : 16061776
          },
          {
            "filename" :
                "/tests/UnitTests/HBLCommands/ExclamationPointEquals.bf",
            "start" : 16061776,
            "end" : 16063987
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/ExecuteAFile.bf",
            "start" : 16063987,
            "end" : 16065430
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/ExecuteCommands.bf",
            "start" : 16065430,
            "end" : 16067089
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Exp.bf",
            "start" : 16067089,
            "end" : 16069570
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Export.bf",
            "start" : 16069570,
            "end" : 16070933
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Ffunction.bf",
            "start" : 16070933,
            "end" : 16071878
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/FindRoot.bf",
            "start" : 16071878,
            "end" : 16074039
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/For.bf",
            "start" : 16074039,
            "end" : 16075496
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Format.bf",
            "start" : 16075496,
            "end" : 16079703
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Function.bf",
            "start" : 16079703,
            "end" : 16081192
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Gamma.bf",
            "start" : 16081192,
            "end" : 16083521
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/GammaDist.bf",
            "start" : 16083521,
            "end" : 16090668
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/GetDataInfo.bf",
            "start" : 16090668,
            "end" : 16092179
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/GetInformation.bf",
            "start" : 16092179,
            "end" : 16096825
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/GetString.bf",
            "start" : 16096825,
            "end" : 16108263
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Greaterthan.bf",
            "start" : 16108263,
            "end" : 16111659
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Greaterthanorequalto.bf",
            "start" : 16111659,
            "end" : 16114742
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/HarvestFrequencies.bf",
            "start" : 16114742,
            "end" : 16118402
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/IBeta.bf",
            "start" : 16118402,
            "end" : 16121958
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/IGamma.bf",
            "start" : 16121958,
            "end" : 16124474
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/If.bf",
            "start" : 16124474,
            "end" : 16125670
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Integrate.bf",
            "start" : 16125670,
            "end" : 16128491
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Inverse.bf",
            "start" : 16128491,
            "end" : 16129161
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/JSON.bf",
            "start" : 16129161,
            "end" : 16130164
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Join.bf",
            "start" : 16130164,
            "end" : 16130669
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/LUDecompose.bf",
            "start" : 16130669,
            "end" : 16132927
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/LUSolve.bf",
            "start" : 16132927,
            "end" : 16135106
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Lessthan.bf",
            "start" : 16135106,
            "end" : 16138468
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Lessthanorequalto.bf",
            "start" : 16138468,
            "end" : 16141409
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/LnGamma.bf",
            "start" : 16141409,
            "end" : 16142882
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/LoadFunctionLibrary.bf",
            "start" : 16142882,
            "end" : 16144430
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Log.bf",
            "start" : 16144430,
            "end" : 16146650
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/MAccess_bracket.bf",
            "start" : 16146650,
            "end" : 16150765
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Max.bf",
            "start" : 16150765,
            "end" : 16153346
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Min.bf",
            "start" : 16153346,
            "end" : 16157225
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Model.bf",
            "start" : 16157225,
            "end" : 16159538
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Percentsign.bf",
            "start" : 16159538,
            "end" : 16162615
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Polynomial.bf",
            "start" : 16162615,
            "end" : 16166354
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Random.bf",
            "start" : 16166354,
            "end" : 16171330
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/ReplicateConstraint.bf",
            "start" : 16171330,
            "end" : 16177252
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/RequireVersion.bf",
            "start" : 16177252,
            "end" : 16177849
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/RerootTree.bf",
            "start" : 16177849,
            "end" : 16181335
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Return.bf",
            "start" : 16181335,
            "end" : 16182735
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Rows.bf",
            "start" : 16182735,
            "end" : 16185030
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Simplex.bf",
            "start" : 16185030,
            "end" : 16185898
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Simplify.bf",
            "start" : 16185898,
            "end" : 16188671
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/SimulateDataSet.bf",
            "start" : 16188671,
            "end" : 16191805
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Sin.bf",
            "start" : 16191805,
            "end" : 16194119
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Sqrt.bf",
            "start" : 16194119,
            "end" : 16196326
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Sscanf.bf",
            "start" : 16196326,
            "end" : 16198738
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Subtract.bf",
            "start" : 16198738,
            "end" : 16200972
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Tan.bf",
            "start" : 16200972,
            "end" : 16203260
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/TestTools.ibf",
            "start" : 16203260,
            "end" : 16204584
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Time.bf",
            "start" : 16204584,
            "end" : 16207585
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/TipCount.bf",
            "start" : 16207585,
            "end" : 16210813
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/TipName.bf",
            "start" : 16210813,
            "end" : 16214191
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Topology.bf",
            "start" : 16214191,
            "end" : 16217286
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Transpose.bf",
            "start" : 16217286,
            "end" : 16219026
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Tree.bf",
            "start" : 16219026,
            "end" : 16221675
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/Type.bf",
            "start" : 16221675,
            "end" : 16222787
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/While.bf",
            "start" : 16222787,
            "end" : 16224167
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/ZCDF.bf",
            "start" : 16224167,
            "end" : 16226442
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/assert.bf",
            "start" : 16226442,
            "end" : 16229009
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/break.bf",
            "start" : 16229009,
            "end" : 16230896
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/fscanf_fprintf.bf",
            "start" : 16230896,
            "end" : 16233412
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/include.bf",
            "start" : 16233412,
            "end" : 16234207
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/lfunction.bf",
            "start" : 16234207,
            "end" : 16235513
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/libv3_iofunctions.bf",
            "start" : 16235513,
            "end" : 16236480
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/libv3_math.bf",
            "start" : 16236480,
            "end" : 16237379
          },
          {
            "filename" :
                "/tests/UnitTests/HBLCommands/libv3_utilityfunctions.bf",
            "start" : 16237379,
            "end" : 16239150
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/namespace.bf",
            "start" : 16239150,
            "end" : 16240856
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/nested/baz.bf",
            "start" : 16240856,
            "end" : 16240877
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/nested/foo.bf",
            "start" : 16240877,
            "end" : 16240890
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/nested/nested2/bar.bf",
            "start" : 16240890,
            "end" : 16240933
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/res/EU3031.nwk",
            "start" : 16240933,
            "end" : 16406112
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/res/SCFG/SCFG.ibf",
            "start" : 16406112,
            "end" : 16418329
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/res/SCFG/scfgG6c.bf",
            "start" : 16418329,
            "end" : 16420920
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/res/SCFG/small.txt",
            "start" : 16420920,
            "end" : 16422373
          },
          {
            "filename" :
                "/tests/UnitTests/HBLCommands/res/replicate_constraint.nex",
            "start" : 16422373,
            "end" : 16458172
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/res/test_likefunc.nex",
            "start" : 16458172,
            "end" : 16493413
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/res/test_likefunc2.nex",
            "start" : 16493413,
            "end" : 16673809
          },
          {
            "filename" :
                "/tests/UnitTests/HBLCommands/testdata/Chinook_Sqlite.sqlite",
            "start" : 16673809,
            "end" : 17765393
          },
          {
            "filename" : "/tests/UnitTests/HBLCommands/tmp/GetURL.txt",
            "start" : 17765393,
            "end" : 17765393
          },
          {
            "filename" : "/tests/data/.DS_Store",
            "start" : 17765393,
            "end" : 17771541
          },
          {
            "filename" : "/tests/data/2.fas",
            "start" : 17771541,
            "end" : 17773295
          },
          {
            "filename" : "/tests/data/2.prot",
            "start" : 17773295,
            "end" : 17773885
          },
          {
            "filename" : "/tests/data/2.tree",
            "start" : 17773885,
            "end" : 17773890
          },
          {
            "filename" : "/tests/data/5.fas",
            "start" : 17773890,
            "end" : 17774265
          },
          {
            "filename" : "/tests/data/CD2-N.nex",
            "start" : 17774265,
            "end" : 17782588
          },
          {
            "filename" : "/tests/data/CD2-tree-with-crazy-lengths.nwk",
            "start" : 17782588,
            "end" : 17782770
          },
          {
            "filename" : "/tests/data/CD2-tree-with-lengths.nwk",
            "start" : 17782770,
            "end" : 17783016
          },
          {
            "filename" : "/tests/data/CD2-tree-without-lengths.nwk",
            "start" : 17783016,
            "end" : 17783134
          },
          {
            "filename" : "/tests/data/CD2.newick",
            "start" : 17783134,
            "end" : 17783380
          },
          {
            "filename" : "/tests/data/CD2.nex",
            "start" : 17783380,
            "end" : 17791703
          },
          {
            "filename" : "/tests/data/CD2.nex.fits",
            "start" : 17791703,
            "end" : 17827507
          },
          {
            "filename" : "/tests/data/CD2.phylip",
            "start" : 17827507,
            "end" : 17835026
          },
          {
            "filename" : "/tests/data/CD2_AA.fna",
            "start" : 17835026,
            "end" : 17837224
          },
          {
            "filename" : "/tests/data/CD2_noTree.nex",
            "start" : 17837224,
            "end" : 17843220
          },
          {
            "filename" : "/tests/data/CD2_reduced.fasta",
            "start" : 17843220,
            "end" : 17843476
          },
          {
            "filename" : "/tests/data/CD2_reduced.fna",
            "start" : 17843476,
            "end" : 17843947
          },
          {
            "filename" : "/tests/data/CD2_reduced.nhx",
            "start" : 17843947,
            "end" : 17844162
          },
          {
            "filename" : "/tests/data/ExecuteAFileSetXTo5.bf",
            "start" : 17844162,
            "end" : 17844168
          },
          {
            "filename" : "/tests/data/HIVenvSweden.seq",
            "start" : 17844168,
            "end" : 17848249
          },
          {
            "filename" : "/tests/data/HMM2_synthetic.fas",
            "start" : 17848249,
            "end" : 17859562
          },
          {
            "filename" : "/tests/data/HMM4_synthetic.fas",
            "start" : 17859562,
            "end" : 17870875
          },
          {
            "filename" : "/tests/data/fluHA.nex",
            "start" : 17870875,
            "end" : 18228132
          },
          {
            "filename" : "/tests/data/fluHA_codon.nex",
            "start" : 18228132,
            "end" : 18590376
          },
          {
            "filename" : "/tests/data/integrase_BDA.nex",
            "start" : 18590376,
            "end" : 18600981
          },
          {
            "filename" : "/tests/data/mammalian_REM2_codons.SA-filtered.fas",
            "start" : 18600981,
            "end" : 18953134
          },
          {
            "filename" : "/tests/data/mammalian_REM2_codons.SA.fasta",
            "start" : 18953134,
            "end" : 19314387
          },
          {
            "filename" : "/tests/data/mtDNA.fas",
            "start" : 19314387,
            "end" : 19991980
          },
          {
            "filename" : "/tests/libv3/.DS_Store",
            "start" : 19991980,
            "end" : 20002224
          },
          {
            "filename" : "/tests/libv3/ABSREL-MH.wbf",
            "start" : 20002224,
            "end" : 20003873
          },
          {
            "filename" : "/tests/libv3/ABSREL.wbf",
            "start" : 20003873,
            "end" : 20005599
          },
          {
            "filename" : "/tests/libv3/BGM.wbf",
            "start" : 20005599,
            "end" : 20007690
          },
          {
            "filename" : "/tests/libv3/BUSTED-MH.wbf",
            "start" : 20007690,
            "end" : 20008530
          },
          {
            "filename" : "/tests/libv3/BUSTED-SRV.wbf",
            "start" : 20008530,
            "end" : 20009381
          },
          {
            "filename" : "/tests/libv3/BUSTED.wbf",
            "start" : 20009381,
            "end" : 20010414
          },
          {
            "filename" : "/tests/libv3/CFEL.wbf",
            "start" : 20010414,
            "end" : 20016924
          },
          {
            "filename" : "/tests/libv3/FADE.wbf",
            "start" : 20016924,
            "end" : 20018636
          },
          {
            "filename" : "/tests/libv3/FEL.wbf",
            "start" : 20018636,
            "end" : 20021176
          },
          {
            "filename" : "/tests/libv3/FMM.wbf",
            "start" : 20021176,
            "end" : 20022904
          },
          {
            "filename" : "/tests/libv3/FUBAR.wbf",
            "start" : 20022904,
            "end" : 20026408
          },
          {
            "filename" : "/tests/libv3/GARD.wbf",
            "start" : 20026408,
            "end" : 20027076
          },
          {
            "filename" : "/tests/libv3/LEISR.wbf",
            "start" : 20027076,
            "end" : 20028094
          },
          {
            "filename" : "/tests/libv3/MEME-partitioned.wbf",
            "start" : 20028094,
            "end" : 20030641
          },
          {
            "filename" : "/tests/libv3/MEME.wbf",
            "start" : 20030641,
            "end" : 20036607
          },
          {
            "filename" : "/tests/libv3/MG94-REV-Constrain-Branch-Lengths.bf",
            "start" : 20036607,
            "end" : 20040572
          },
          {
            "filename" : "/tests/libv3/ProteinGTRFit.wbf",
            "start" : 20040572,
            "end" : 20042024
          },
          {
            "filename" : "/tests/libv3/RELAX.wbf",
            "start" : 20042024,
            "end" : 20043996
          },
          {
            "filename" : "/tests/libv3/SLAC-partitioned.wbf",
            "start" : 20043996,
            "end" : 20045115
          },
          {
            "filename" : "/tests/libv3/SLAC.wbf",
            "start" : 20045115,
            "end" : 20046223
          },
          {
            "filename" : "/tests/libv3/algae-mtDNA.wbf",
            "start" : 20046223,
            "end" : 20047134
          },
          {
            "filename" : "/tests/libv3/alignments.bf",
            "start" : 20047134,
            "end" : 20047729
          },
          {
            "filename" : "/tests/libv3/ciliate-code.wbf",
            "start" : 20047729,
            "end" : 20048606
          },
          {
            "filename" : "/tests/libv3/data/.DS_Store",
            "start" : 20048606,
            "end" : 20054754
          },
          {
            "filename" : "/tests/libv3/data/.gitignore",
            "start" : 20054754,
            "end" : 20054786
          },
          {
            "filename" : "/tests/libv3/data/CD2-filtered.nex",
            "start" : 20054786,
            "end" : 20056646
          },
          {
            "filename" : "/tests/libv3/data/CD2.fits",
            "start" : 20056646,
            "end" : 20080389
          },
          {
            "filename" : "/tests/libv3/data/CD2.nex",
            "start" : 20080389,
            "end" : 20088544
          },
          {
            "filename" : "/tests/libv3/data/CD2.prot",
            "start" : 20088544,
            "end" : 20090952
          },
          {
            "filename" : "/tests/libv3/data/CD2_AA.fna",
            "start" : 20090952,
            "end" : 20093150
          },
          {
            "filename" : "/tests/libv3/data/COXI.nex",
            "start" : 20093150,
            "end" : 20126589
          },
          {
            "filename" : "/tests/libv3/data/ENSG00000001629.fa",
            "start" : 20126589,
            "end" : 20206437
          },
          {
            "filename" : "/tests/libv3/data/HRVI.nex",
            "start" : 20206437,
            "end" : 20212085
          },
          {
            "filename" : "/tests/libv3/data/adh.nex",
            "start" : 20212085,
            "end" : 20230377
          },
          {
            "filename" : "/tests/libv3/data/algae_coxI_aligned.fas",
            "start" : 20230377,
            "end" : 20243691
          },
          {
            "filename" : "/tests/libv3/data/ciliates.fas",
            "start" : 20243691,
            "end" : 20253243
          },
          {
            "filename" : "/tests/libv3/data/dat.nuc",
            "start" : 20253243,
            "end" : 20253921
          },
          {
            "filename" : "/tests/libv3/data/dat.prot",
            "start" : 20253921,
            "end" : 20254355
          },
          {
            "filename" : "/tests/libv3/data/json-read.py",
            "start" : 20254355,
            "end" : 20254704
          },
          {
            "filename" : "/tests/libv3/data/mammals.tre",
            "start" : 20254704,
            "end" : 20255369
          },
          {
            "filename" : "/tests/libv3/data/ncov.fasta",
            "start" : 20255369,
            "end" : 20275158
          },
          {
            "filename" : "/tests/libv3/data/partitioned.nex",
            "start" : 20275158,
            "end" : 20312318
          },
          {
            "filename" :
                "/tests/libv3/data/partitioned.nex.EtaqIXyn.likelihoodFunction.bf",
            "start" : 20312318,
            "end" : 20415252
          },
          {
            "filename" : "/tests/libv3/data/partitioned.nex.SLAC.bf",
            "start" : 20415252,
            "end" : 20518003
          },
          {
            "filename" :
                "/tests/libv3/data/protgtr_fitter_alignments/prot1.dat",
            "start" : 20518003,
            "end" : 20528822
          },
          {
            "filename" :
                "/tests/libv3/data/protgtr_fitter_alignments/prot2.dat",
            "start" : 20528822,
            "end" : 20532891
          },
          {
            "filename" :
                "/tests/libv3/data/protgtr_fitter_alignments/prot3.dat",
            "start" : 20532891,
            "end" : 20541679
          },
          {
            "filename" :
                "/tests/libv3/data/protgtr_fitter_alignments/prot4.dat",
            "start" : 20541679,
            "end" : 20547585
          },
          {
            "filename" :
                "/tests/libv3/data/protgtr_fitter_alignments/prot5.dat",
            "start" : 20547585,
            "end" : 20560527
          },
          {
            "filename" : "/tests/libv3/data/protgtr_fitter_lines_raw.txt",
            "start" : 20560527,
            "end" : 20560836
          },
          {
            "filename" : "/tests/libv3/models-codon.bf",
            "start" : 20560836,
            "end" : 20562098
          },
          {
            "filename" : "/tests/libv3/models/absrel.bf",
            "start" : 20562098,
            "end" : 20568353
          },
          {
            "filename" : "/tests/libv3/models/hky85.bf",
            "start" : 20568353,
            "end" : 20571286
          },
          {
            "filename" : "/tests/libv3/mtDNA-code.wbf",
            "start" : 20571286,
            "end" : 20572175
          },
          {
            "filename" : "/tests/libv3/rate-variation.bf",
            "start" : 20572175,
            "end" : 20577780
          },
          {
            "filename" : "/tests/libv3/relative_nucleotide_rates.wbf",
            "start" : 20577780,
            "end" : 20577914
          },
          {
            "filename" : "/tests/libv3/relative_prot_rates.wbf",
            "start" : 20577914,
            "end" : 20578091
          },
          {
            "filename" : "/tests/libv3/shared.bf",
            "start" : 20578091,
            "end" : 20578315
          },
          {
            "filename" : "/tests/libv3/support/FitMG94.bf",
            "start" : 20578315,
            "end" : 20590879
          },
          {
            "filename" : "/tests/libv3/trees.bf",
            "start" : 20590879,
            "end" : 20596084
          }
        ],
        "remote_package_size" : 20596084
      });
    })();

    // end include:
    // /var/folders/_6/5pbbnj3j5x1b_vnpy1q7fzlc0000gp/T/tmphcfvyuxw.js include:
    // /var/folders/_6/5pbbnj3j5x1b_vnpy1q7fzlc0000gp/T/tmpq0e2gmpv.js

    // All the pre-js content up to here must remain later on, we need to run
    // it.
    if ((typeof ENVIRONMENT_IS_WASM_WORKER != 'undefined' &&
         ENVIRONMENT_IS_WASM_WORKER) ||
        (typeof ENVIRONMENT_IS_PTHREAD != 'undefined' &&
         ENVIRONMENT_IS_PTHREAD) ||
        (typeof ENVIRONMENT_IS_AUDIO_WORKLET != 'undefined' &&
         ENVIRONMENT_IS_AUDIO_WORKLET))
      Module['preRun'] = [];
    var necessaryPreJSTasks = Module['preRun'].slice();
    // end include:
    // /var/folders/_6/5pbbnj3j5x1b_vnpy1q7fzlc0000gp/T/tmpq0e2gmpv.js
    // include: /var/folders/_6/5pbbnj3j5x1b_vnpy1q7fzlc0000gp/T/tmpc40rod77.js

    if (!Module['preRun'])
      throw 'Module.preRun should exist because file support used it; did a pre-js delete it?';
    necessaryPreJSTasks.forEach((task) => {
      if (Module['preRun'].indexOf(task) < 0)
        throw 'All preRun tasks that exist before user pre-js code should remain after; did you replace Module or modify Module.preRun?';
    });
    // end include:
    // /var/folders/_6/5pbbnj3j5x1b_vnpy1q7fzlc0000gp/T/tmpc40rod77.js

    var arguments_ = [];
    var thisProgram = './this.program';
    var quit_ = (status, toThrow) => { throw toThrow; };

    if (ENVIRONMENT_IS_WORKER) {
      _scriptName = self.location.href;
    }

    // `/` should be present at the end if `scriptDirectory` is not empty
    var scriptDirectory = '';
    function locateFile(path) {
      if (Module['locateFile']) {
        return Module['locateFile'](path, scriptDirectory);
      }
      return scriptDirectory + path;
    }

    // Hooks that are implemented differently in different runtime environments.
    var readAsync, readBinary;

    if (ENVIRONMENT_IS_SHELL) {

      const isNode = typeof process == 'object' && process.versions?.node &&
                     process.type != 'renderer';
      if (isNode || typeof window == 'object' ||
          typeof WorkerGlobalScope != 'undefined')
        throw new Error(
            'not compiled for this environment (did you build to HTML and try to run it not on the web, or set ENVIRONMENT to something - like node - and run it someplace else - like on the web?)');

    } else

      // Note that this includes Node.js workers when relevant (pthreads is
      // enabled). Node.js workers are detected as a combination of
      // ENVIRONMENT_IS_WORKER and ENVIRONMENT_IS_NODE.
      if (ENVIRONMENT_IS_WEB || ENVIRONMENT_IS_WORKER) {
        try {
          scriptDirectory =
              new URL('.', _scriptName).href; // includes trailing slash
        } catch {
          // Must be a `blob:` or `data:` URL (e.g.
          // `blob:http://site.com/etc/etc`), we cannot infer anything from
          // them.
        }

        if (!(typeof window == 'object' ||
              typeof WorkerGlobalScope != 'undefined'))
          throw new Error(
              'not compiled for this environment (did you build to HTML and try to run it not on the web, or set ENVIRONMENT to something - like node - and run it someplace else - like on the web?)');

        {
          // include: web_or_worker_shell_read.js
          if (ENVIRONMENT_IS_WORKER) {
            readBinary = (url) => {
              var xhr = new XMLHttpRequest();
              xhr.open('GET', url, false);
              xhr.responseType = 'arraybuffer';
              xhr.send(null);
              return new Uint8Array(/** @type{!ArrayBuffer} */ (xhr.response));
            };
          }

          readAsync = async (url) => {
            assert(!isFileURI(url),
                   "readAsync does not work with file:// URLs");
            var response = await fetch(url, {credentials : 'same-origin'});
            if (response.ok) {
              return response.arrayBuffer();
            }
            throw new Error(response.status + ' : ' + response.url);
          };
          // end include: web_or_worker_shell_read.js
        }
      } else {
        throw new Error('environment detection error');
      }

    var out = console.log.bind(console);
    var err = console.error.bind(console);

    var IDBFS = 'IDBFS is no longer included by default; build with -lidbfs.js';

    var FETCHFS =
        'FETCHFS is no longer included by default; build with -lfetchfs.js';
    var ICASEFS =
        'ICASEFS is no longer included by default; build with -licasefs.js';
    var JSFILEFS =
        'JSFILEFS is no longer included by default; build with -ljsfilefs.js';
    var OPFS = 'OPFS is no longer included by default; build with -lopfs.js';

    var NODEFS =
        'NODEFS is no longer included by default; build with -lnodefs.js';

    // perform assertions in shell.js after we set up out() and err(), as
    // otherwise if an assertion fails it cannot print the message

    assert(
        !ENVIRONMENT_IS_NODE,
        'node environment detected but not enabled at build time.  Add `node` to `-sENVIRONMENT` to enable.');

    assert(
        !ENVIRONMENT_IS_SHELL,
        'shell environment detected but not enabled at build time.  Add `shell` to `-sENVIRONMENT` to enable.');

    // end include: shell.js

    // include: preamble.js
    // === Preamble library stuff ===

    // Documentation for the public APIs defined in this file must be updated
    // in:
    //    site/source/docs/api_reference/preamble.js.rst
    // A prebuilt local version of the documentation is available at:
    //    site/build/text/docs/api_reference/preamble.js.txt
    // You can also build docs locally as HTML or other formats in site/
    // An online HTML version (which may be of a different version of
    // Emscripten)
    //    is up at
    //    http://kripken.github.io/emscripten-site/docs/api_reference/preamble.js.html

    var wasmBinary;

    if (typeof WebAssembly != 'object') {
      err('no native wasm support detected');
    }

    // Wasm globals

    //========================================
    // Runtime essentials
    //========================================

    // whether we are quitting the application. no code should run after this.
    // set in exit() and abort()
    var ABORT = false;

    // set by exit() and abort().  Passed to 'onExit' handler.
    // NOTE: This is also used as the process return code code in shell
    // environments but only when noExitRuntime is false.
    var EXITSTATUS;

    // In STRICT mode, we only define assert() when ASSERTIONS is set.  i.e. we
    // don't define it at all in release modes.  This matches the behaviour of
    // MINIMAL_RUNTIME.
    // TODO(sbc): Make this the default even without STRICT enabled.
    /** @type {function(*, string=)} */
    function assert(condition, text) {
      if (!condition) {
        abort('Assertion failed' + (text ? ': ' + text : ''));
      }
    }

    // We used to include malloc/free by default in the past. Show a helpful
    // error in builds with assertions.
    function _malloc() {
      abort(
          'malloc() called but not included in the build - add `_malloc` to EXPORTED_FUNCTIONS');
    }

    /**
     * Indicates whether filename is delivered via file protocol (as opposed to
     * http/https)
     * @noinline
     */
    var isFileURI = (filename) => filename.startsWith('file://');

    // include: runtime_common.js
    // include: runtime_stack_check.js
    // Initializes the stack cookie. Called at the startup of main and at the
    // startup of each thread in pthreads mode.
    function writeStackCookie() {
      var max = _emscripten_stack_get_end();
      assert((max & 3) == 0);
      // If the stack ends at address zero we write our cookies 4 bytes into the
      // stack.  This prevents interference with SAFE_HEAP and ASAN which also
      // monitor writes to address zero.
      if (max == 0) {
        max += 4;
      }
      // The stack grow downwards towards _emscripten_stack_get_end.
      // We write cookies to the final two words in the stack and detect if they
      // are ever overwritten.
      HEAPU32[((max) >> 2)] = 0x02135467;
      HEAPU32[(((max) + (4)) >> 2)] = 0x89BACDFE;
      // Also test the global address 0 for integrity.
      HEAPU32[((0) >> 2)] = 1668509029;
    }

    function checkStackCookie() {
      if (ABORT)
        return;
      var max = _emscripten_stack_get_end();
      // See writeStackCookie().
      if (max == 0) {
        max += 4;
      }
      var cookie1 = HEAPU32[((max) >> 2)];
      var cookie2 = HEAPU32[(((max) + (4)) >> 2)];
      if (cookie1 != 0x02135467 || cookie2 != 0x89BACDFE) {
        abort(`Stack overflow! Stack cookie has been overwritten at ${
            ptrToString(
                max)}, expected hex dwords 0x89BACDFE and 0x2135467, but received ${
            ptrToString(cookie2)} ${ptrToString(cookie1)}`);
      }
      // Also test the global address 0 for integrity.
      if (HEAPU32[((0) >> 2)] != 0x63736d65 /* 'emsc' */) {
        abort(
            'Runtime error: The application has corrupted its heap memory area (address zero)!');
      }
    }
    // end include: runtime_stack_check.js
    // include: runtime_exceptions.js
    // end include: runtime_exceptions.js
    // include: runtime_debug.js
    var runtimeDebug = true; // Switch to false at runtime to disable logging at
                             // the right times

    // Used by XXXXX_DEBUG settings to output debug messages.
    function dbg(...args) {
      if (!runtimeDebug && typeof runtimeDebug != 'undefined')
        return;
      // TODO(sbc): Make this configurable somehow.  Its not always convenient
      // for logging to show up as warnings.
      console.warn(...args);
    }

    // Endianness check
    (() => {
      var h16 = new Int16Array(1);
      var h8 = new Int8Array(h16.buffer);
      h16[0] = 0x6373;
      if (h8[0] !== 0x73 || h8[1] !== 0x63)
        throw 'Runtime error: expected the system to be little-endian! (Run with -sSUPPORT_BIG_ENDIAN to bypass)';
    })();

    function consumedModuleProp(prop) {
      if (!Object.getOwnPropertyDescriptor(Module, prop)) {
        Object.defineProperty(Module, prop, {
          configurable : true,
          set() {
            abort(`Attempt to set \`Module.${
                prop}\` after it has already been processed.  This can happen, for example, when code is injected via '--post-js' rather than '--pre-js'`);
          }
        });
      }
    }

    function makeInvalidEarlyAccess(name) {
      return () => assert(
                 false,
                 `call to '${
                     name}' via reference taken before Wasm module initialization`);
    }

    function ignoredModuleProp(prop) {
      if (Object.getOwnPropertyDescriptor(Module, prop)) {
        abort(`\`Module.${prop}\` was supplied but \`${
            prop}\` not included in INCOMING_MODULE_JS_API`);
      }
    }

    // forcing the filesystem exports a few things by default
    function isExportedByForceFilesystem(name) {
      return name === 'FS_createPath' || name === 'FS_createDataFile' ||
             name === 'FS_createPreloadedFile' || name === 'FS_preloadFile' ||
             name === 'FS_unlink' || name === 'addRunDependency' ||
             // The old FS has some functionality that WasmFS lacks.
             name === 'FS_createLazyFile' || name === 'FS_createDevice' ||
             name === 'removeRunDependency';
    }

    /**
     * Intercept access to a global symbol.  This enables us to give informative
     * warnings/errors when folks attempt to use symbols they did not include in
     * their build, or no symbols that no longer exist.
     */
    function hookGlobalSymbolAccess(sym, func) {
      // In MODULARIZE mode the generated code runs inside a function scope and
      // not the global scope, and JavaScript does not provide access to
      // function scopes so we cannot dynamically modify the scrope using
      // `defineProperty` in this case.
      //
      // In this mode we simply ignore requests for `hookGlobalSymbolAccess`.
      // Since this is a debug-only feature, skipping it is not major issue.
    }

    function missingGlobal(sym, msg) {
      hookGlobalSymbolAccess(sym, () => {
        warnOnce(`\`${sym}\` is not longer defined by emscripten. ${msg}`);
      });
    }

    missingGlobal('buffer', 'Please use HEAP8.buffer or wasmMemory.buffer');
    missingGlobal('asm', 'Please use wasmExports instead');

    function missingLibrarySymbol(sym) {
      hookGlobalSymbolAccess(sym, () => {
        // Can't `abort()` here because it would break code that does runtime
        // checks.  e.g. `if (typeof SDL === 'undefined')`.
        var msg = `\`${
            sym}\` is a library symbol and not included by default; add it to your library.js __deps or to DEFAULT_LIBRARY_FUNCS_TO_INCLUDE on the command line`;
        // DEFAULT_LIBRARY_FUNCS_TO_INCLUDE requires the name as it appears in
        // library.js, which means $name for a JS name with no prefix, or name
        // for a JS name like _name.
        var librarySymbol = sym;
        if (!librarySymbol.startsWith('_')) {
          librarySymbol = '$' + sym;
        }
        msg += ` (e.g. -sDEFAULT_LIBRARY_FUNCS_TO_INCLUDE='${librarySymbol}')`;
        if (isExportedByForceFilesystem(sym)) {
          msg +=
              '. Alternatively, forcing filesystem support (-sFORCE_FILESYSTEM) can export this for you';
        }
        warnOnce(msg);
      });

      // Any symbol that is not included from the JS library is also (by
      // definition) not exported on the Module object.
      unexportedRuntimeSymbol(sym);
    }

    function unexportedRuntimeSymbol(sym) {
      if (!Object.getOwnPropertyDescriptor(Module, sym)) {
        Object.defineProperty(Module, sym, {
          configurable : true,
          get() {
            var msg = `'${
                sym}' was not exported. add it to EXPORTED_RUNTIME_METHODS (see the Emscripten FAQ)`;
            if (isExportedByForceFilesystem(sym)) {
              msg +=
                  '. Alternatively, forcing filesystem support (-sFORCE_FILESYSTEM) can export this for you';
            }
            abort(msg);
          }
        });
      }
    }

    // end include: runtime_debug.js
    var readyPromiseResolve, readyPromiseReject;

    // Memory management

    var wasmMemory;

    var
        /** @type {!Int8Array} */
        HEAP8,
        /** @type {!Uint8Array} */
        HEAPU8,
        /** @type {!Int16Array} */
        HEAP16,
        /** @type {!Uint16Array} */
        HEAPU16,
        /** @type {!Int32Array} */
        HEAP32,
        /** @type {!Uint32Array} */
        HEAPU32,
        /** @type {!Float32Array} */
        HEAPF32,
        /** @type {!Float64Array} */
        HEAPF64;

    // BigInt64Array type is not correctly defined in closure
    var
        /** not-@type {!BigInt64Array} */
        HEAP64,
        /* BigUint64Array type is not correctly defined in closure
        /** not-@type {!BigUint64Array} */
        HEAPU64;

    var runtimeInitialized = false;

    function updateMemoryViews() {
      var b = wasmMemory.buffer;
      HEAP8 = new Int8Array(b);
      HEAP16 = new Int16Array(b);
      HEAPU8 = new Uint8Array(b);
      HEAPU16 = new Uint16Array(b);
      HEAP32 = new Int32Array(b);
      HEAPU32 = new Uint32Array(b);
      HEAPF32 = new Float32Array(b);
      HEAPF64 = new Float64Array(b);
      HEAP64 = new BigInt64Array(b);
      HEAPU64 = new BigUint64Array(b);
    }

    // include: memoryprofiler.js
    // end include: memoryprofiler.js
    // end include: runtime_common.js
    assert(typeof Int32Array != 'undefined' &&
               typeof Float64Array !== 'undefined' &&
               Int32Array.prototype.subarray != undefined &&
               Int32Array.prototype.set != undefined,
           'JS engine does not provide full typed array support');

    function preRun() {
      if (Module['preRun']) {
        if (typeof Module['preRun'] == 'function')
          Module['preRun'] = [ Module['preRun'] ];
        while (Module['preRun'].length) {
          addOnPreRun(Module['preRun'].shift());
        }
      }
      consumedModuleProp('preRun');
      // Begin ATPRERUNS hooks
      callRuntimeCallbacks(onPreRuns);
      // End ATPRERUNS hooks
    }

    function initRuntime() {
      assert(!runtimeInitialized);
      runtimeInitialized = true;

      checkStackCookie();

      // Begin ATINITS hooks
      if (!Module['noFSInit'] && !FS.initialized)
        FS.init();
      TTY.init();
      // End ATINITS hooks

      wasmExports['__wasm_call_ctors']();

      // Begin ATPOSTCTORS hooks
      FS.ignorePermissions = false;
      // End ATPOSTCTORS hooks
    }

    function preMain() {
      checkStackCookie();
      // No ATMAINS hooks
    }

    function postRun() {
      checkStackCookie();
      // PThreads reuse the runtime from the main thread.

      if (Module['postRun']) {
        if (typeof Module['postRun'] == 'function')
          Module['postRun'] = [ Module['postRun'] ];
        while (Module['postRun'].length) {
          addOnPostRun(Module['postRun'].shift());
        }
      }
      consumedModuleProp('postRun');

      // Begin ATPOSTRUNS hooks
      callRuntimeCallbacks(onPostRuns);
      // End ATPOSTRUNS hooks
    }

    // A counter of dependencies for calling run(). If we need to
    // do asynchronous work before running, increment this and
    // decrement it. Incrementing must happen in a place like
    // Module.preRun (used by emcc to add file preloading).
    // Note that you can add dependencies in preRun, even though
    // it happens right before run - run will be postponed until
    // the dependencies are met.
    var runDependencies = 0;
    var dependenciesFulfilled = null; // overridden to take different actions
                                      // when all run dependencies are fulfilled
    var runDependencyTracking = {};
    var runDependencyWatcher = null;

    function addRunDependency(id) {
      runDependencies++;

      Module['monitorRunDependencies']?.(runDependencies);

      assert(id, 'addRunDependency requires an ID')
      assert(!runDependencyTracking[id]);
      runDependencyTracking[id] = 1;
      if (runDependencyWatcher === null && typeof setInterval != 'undefined') {
        // Check for missing dependencies every few seconds
        runDependencyWatcher = setInterval(() => {
          if (ABORT) {
            clearInterval(runDependencyWatcher);
            runDependencyWatcher = null;
            return;
          }
          var shown = false;
          for (var dep in runDependencyTracking) {
            if (!shown) {
              shown = true;
              err('still waiting on run dependencies:');
            }
            err(`dependency: ${dep}`);
          }
          if (shown) {
            err('(end of list)');
          }
        }, 10000);
      }
    }

    function removeRunDependency(id) {
      runDependencies--;

      Module['monitorRunDependencies']?.(runDependencies);

      assert(id, 'removeRunDependency requires an ID');
      assert(runDependencyTracking[id]);
      delete runDependencyTracking[id];
      if (runDependencies == 0) {
        if (runDependencyWatcher !== null) {
          clearInterval(runDependencyWatcher);
          runDependencyWatcher = null;
        }
        if (dependenciesFulfilled) {
          var callback = dependenciesFulfilled;
          dependenciesFulfilled = null;
          callback(); // can add another dependenciesFulfilled
        }
      }
    }

    /** @param {string|number=} what */
    function abort(what) {
      Module['onAbort']?.(what);

      what = 'Aborted(' + what + ')';
      // TODO(sbc): Should we remove printing and leave it up to whoever
      // catches the exception?
      err(what);

      ABORT = true;

      // Use a wasm runtime error, because a JS error might be seen as a foreign
      // exception, which means we'd run destructors on it. We need the error to
      // simply make the program stop.
      // FIXME This approach does not work in Wasm EH because it currently does
      // not assume all RuntimeErrors are from traps; it decides whether a
      // RuntimeError is from a trap or not based on a hidden field within the
      // object. So at the moment we don't have a way of throwing a wasm trap
      // from JS. TODO Make a JS API that allows this in the wasm spec.

      // Suppress closure compiler warning here. Closure compiler's builtin
      // extern definition for WebAssembly.RuntimeError claims it takes no
      // arguments even though it can.
      // TODO(https://github.com/google/closure-compiler/pull/3913): Remove
      // if/when upstream closure gets fixed. See above, in the meantime, we
      // resort to wasm code for trapping.
      //
      // In case abort() is called before the module is initialized, wasmExports
      // and its exported '__trap' function is not available, in which case we
      // throw a RuntimeError.
      //
      // We trap instead of throwing RuntimeError to prevent infinite-looping in
      // Wasm EH code (because RuntimeError is considered as a foreign exception
      // and caught by 'catch_all'), but in case throwing RuntimeError is fine
      // because the module has not even been instantiated, even less running.
      if (runtimeInitialized) {
        ___trap();
      }
      /** @suppress {checkTypes} */
      var e = new WebAssembly.RuntimeError(what);

      readyPromiseReject?.(e);
      // Throw the error whether or not MODULARIZE is set because abort is used
      // in code paths apart from instantiation where an exception is expected
      // to be thrown when abort is called.
      throw e;
    }

    function createExportWrapper(name, nargs) {
      return (...args) => {
        assert(
            runtimeInitialized,
            `native function \`${name}\` called before runtime initialization`);
        var f = wasmExports[name];
        assert(f, `exported native function \`${name}\` not found`);
        // Only assert for too many arguments. Too few can be valid since the
        // missing arguments will be zero filled.
        assert(args.length <= nargs,
               `native function \`${name}\` called with ${
                   args.length} args but expects ${nargs}`);
        return f(...args);
      };
    }

    var wasmBinaryFile;

    function findWasmBinary() { return locateFile('hyphy.wasm'); }

    function getBinarySync(file) {
      if (file == wasmBinaryFile && wasmBinary) {
        return new Uint8Array(wasmBinary);
      }
      if (readBinary) {
        return readBinary(file);
      }
      throw 'both async and sync fetching of the wasm failed';
    }

    async function getWasmBinary(binaryFile) {
      // If we don't have the binary yet, load it asynchronously using
      // readAsync.
      if (!wasmBinary) {
        // Fetch the binary using readAsync
        try {
          var response = await readAsync(binaryFile);
          return new Uint8Array(response);
        } catch {
          // Fall back to getBinarySync below;
        }
      }

      // Otherwise, getBinarySync should be able to get it synchronously
      return getBinarySync(binaryFile);
    }

    async function instantiateArrayBuffer(binaryFile, imports) {
      try {
        var binary = await getWasmBinary(binaryFile);
        var instance = await WebAssembly.instantiate(binary, imports);
        return instance;
      } catch (reason) {
        err(`failed to asynchronously prepare wasm: ${reason}`);

        // Warn on some common problems.
        if (isFileURI(wasmBinaryFile)) {
          err(`warning: Loading from a file URI (${
              wasmBinaryFile}) is not supported in most browsers. See https://emscripten.org/docs/getting_started/FAQ.html#how-do-i-run-a-local-webserver-for-testing-why-does-my-program-stall-in-downloading-or-preparing`);
        }
        abort(reason);
      }
    }

    async function instantiateAsync(binary, binaryFile, imports) {
      if (!binary) {
        try {
          var response = fetch(binaryFile, {credentials : 'same-origin'});
          var instantiationResult =
              await WebAssembly.instantiateStreaming(response, imports);
          return instantiationResult;
        } catch (reason) {
          // We expect the most common failure cause to be a bad MIME type for
          // the binary, in which case falling back to ArrayBuffer instantiation
          // should work.
          err(`wasm streaming compile failed: ${reason}`);
          err('falling back to ArrayBuffer instantiation');
          // fall back of instantiateArrayBuffer below
        };
      }
      return instantiateArrayBuffer(binaryFile, imports);
    }

    function getWasmImports() {
      // prepare imports
      return { 'env': wasmImports, 'wasi_snapshot_preview1': wasmImports, }
    }

    // Create the wasm instance.
    // Receives the wasm imports, returns the exports.
    async function createWasm() {
      // Load the wasm module and create an instance of using native support in
      // the JS engine. handle a generated wasm instance, receiving its exports
      // and performing other necessary setup
      /** @param {WebAssembly.Module=} module*/
      function receiveInstance(instance, module) {
        wasmExports = instance.exports;

        wasmMemory = wasmExports['memory'];

        assert(wasmMemory, 'memory not found in wasm exports');
        updateMemoryViews();

        wasmTable = wasmExports['__indirect_function_table'];

        assert(wasmTable, 'table not found in wasm exports');

        ___cpp_exception = wasmExports['__cpp_exception'];
        ;

        assignWasmExports(wasmExports);
        removeRunDependency('wasm-instantiate');
        return wasmExports;
      }
      // wait for the pthread pool (if any)
      addRunDependency('wasm-instantiate');

      // Prefer streaming instantiation if available.
      // Async compilation can be confusing when an error on the page overwrites
      // Module (for example, if the order of elements is wrong, and the one
      // defining Module is later), so we save Module and check it later.
      var trueModule = Module;
      function receiveInstantiationResult(result) {
        // 'result' is a ResultObject object which has both the module and
        // instance. receiveInstance() will swap in the exports (to Module.asm)
        // so they can be called
        assert(
            Module === trueModule,
            'the Module object should not be replaced during async compilation - perhaps the order of HTML elements is wrong?');
        trueModule = null;
        // TODO: Due to Closure regression
        // https://github.com/google/closure-compiler/issues/3193, the above
        // line no longer optimizes out down to the following line. When the
        // regression is fixed, can restore the above PTHREADS-enabled path.
        return receiveInstance(result['instance']);
      }

      var info = getWasmImports();

      // User shell pages can write their own Module.instantiateWasm =
      // function(imports, successCallback) callback to manually instantiate the
      // Wasm module themselves. This allows pages to run the instantiation
      // parallel to any other async startup actions they are performing. Also
      // pthreads and wasm workers initialize the wasm instance through this
      // path.
      if (Module['instantiateWasm']) {
        return new Promise((resolve, reject) => {
          try {
            Module['instantiateWasm'](
                info, (mod, inst) => { resolve(receiveInstance(mod, inst)); });
          } catch (e) {
            err(`Module.instantiateWasm callback failed with error: ${e}`);
            reject(e);
          }
        });
      }

      wasmBinaryFile ??= findWasmBinary();
      var result = await instantiateAsync(wasmBinary, wasmBinaryFile, info);
      var exports = receiveInstantiationResult(result);
      return exports;
    }

    // end include: preamble.js

    // Begin JS library code

    class ExitStatus {
      name = 'ExitStatus';
      constructor(status) {
        this.message = `Program terminated with exit(${status})`;
        this.status = status;
      }
    }

    var callRuntimeCallbacks = (callbacks) => {
      while (callbacks.length > 0) {
        // Pass the module as the first argument.
        callbacks.shift()(Module);
      }
    };
    var onPostRuns = [];
    var addOnPostRun = (cb) => onPostRuns.push(cb);

    var onPreRuns = [];
    var addOnPreRun = (cb) => onPreRuns.push(cb);

    /**
     * @param {number} ptr
     * @param {string} type
     */
    function getValue(ptr, type = 'i8') {
      if (type.endsWith('*'))
        type = '*';
      switch (type) {
      case 'i1':
        return HEAP8[ptr];
      case 'i8':
        return HEAP8[ptr];
      case 'i16':
        return HEAP16[((ptr) >> 1)];
      case 'i32':
        return HEAP32[((ptr) >> 2)];
      case 'i64':
        return HEAP64[((ptr) >> 3)];
      case 'float':
        return HEAPF32[((ptr) >> 2)];
      case 'double':
        return HEAPF64[((ptr) >> 3)];
      case '*':
        return HEAPU32[((ptr) >> 2)];
      default:
        abort(`invalid type for getValue: ${type}`);
      }
    }

    var noExitRuntime = true;

    var ptrToString = (ptr) => {
      assert(typeof ptr === 'number');
      // With CAN_ADDRESS_2GB or MEMORY64, pointers are already unsigned.
      ptr >>>= 0;
      return '0x' + ptr.toString(16).padStart(8, '0');
    };

    /**
     * @param {number} ptr
     * @param {number} value
     * @param {string} type
     */
    function setValue(ptr, value, type = 'i8') {
      if (type.endsWith('*'))
        type = '*';
      switch (type) {
      case 'i1':
        HEAP8[ptr] = value;
        break;
      case 'i8':
        HEAP8[ptr] = value;
        break;
      case 'i16':
        HEAP16[((ptr) >> 1)] = value;
        break;
      case 'i32':
        HEAP32[((ptr) >> 2)] = value;
        break;
      case 'i64':
        HEAP64[((ptr) >> 3)] = BigInt(value);
        break;
      case 'float':
        HEAPF32[((ptr) >> 2)] = value;
        break;
      case 'double':
        HEAPF64[((ptr) >> 3)] = value;
        break;
      case '*':
        HEAPU32[((ptr) >> 2)] = value;
        break;
      default:
        abort(`invalid type for setValue: ${type}`);
      }
    }

    var warnOnce = (text) => {
      warnOnce.shown ||= {};
      if (!warnOnce.shown[text]) {
        warnOnce.shown[text] = 1;
        err(text);
      }
    };

    var UTF8Decoder =
        typeof TextDecoder != 'undefined' ? new TextDecoder() : undefined;

    var findStringEnd = (heapOrArray, idx, maxBytesToRead, ignoreNul) => {
      var maxIdx = idx + maxBytesToRead;
      if (ignoreNul)
        return maxIdx;
      // TextDecoder needs to know the byte length in advance, it doesn't stop
      // on null terminator by itself. As a tiny code save trick, compare idx
      // against maxIdx using a negation, so that maxBytesToRead=undefined/NaN
      // means Infinity.
      while (heapOrArray[idx] && !(idx >= maxIdx))
        ++idx;
      return idx;
    };

    /**
     * Given a pointer 'idx' to a null-terminated UTF8-encoded string in the
     * given array that contains uint8 values, returns a copy of that string as
     * a Javascript String object. heapOrArray is either a regular array, or a
     * JavaScript typed array view.
     * @param {number=} idx
     * @param {number=} maxBytesToRead
     * @param {boolean=} ignoreNul - If true, the function will not stop on a
     *     NUL character.
     * @return {string}
     */
    var UTF8ArrayToString = (heapOrArray, idx = 0, maxBytesToRead,
                             ignoreNul) => {
      var endPtr = findStringEnd(heapOrArray, idx, maxBytesToRead, ignoreNul);

      // When using conditional TextDecoder, skip it for short strings as the
      // overhead of the native call is not worth it.
      if (endPtr - idx > 16 && heapOrArray.buffer && UTF8Decoder) {
        return UTF8Decoder.decode(heapOrArray.subarray(idx, endPtr));
      }
      var str = '';
      while (idx < endPtr) {
        // For UTF8 byte structure, see:
        // http://en.wikipedia.org/wiki/UTF-8#Description
        // https://www.ietf.org/rfc/rfc2279.txt
        // https://tools.ietf.org/html/rfc3629
        var u0 = heapOrArray[idx++];
        if (!(u0 & 0x80)) {
          str += String.fromCharCode(u0);
          continue;
        }
        var u1 = heapOrArray[idx++] & 63;
        if ((u0 & 0xE0) == 0xC0) {
          str += String.fromCharCode(((u0 & 31) << 6) | u1);
          continue;
        }
        var u2 = heapOrArray[idx++] & 63;
        if ((u0 & 0xF0) == 0xE0) {
          u0 = ((u0 & 15) << 12) | (u1 << 6) | u2;
        } else {
          if ((u0 & 0xF8) != 0xF0)
            warnOnce(
                'Invalid UTF-8 leading byte ' + ptrToString(u0) +
                ' encountered when deserializing a UTF-8 string in wasm memory to a JS string!');
          u0 = ((u0 & 7) << 18) | (u1 << 12) | (u2 << 6) |
               (heapOrArray[idx++] & 63);
        }

        if (u0 < 0x10000) {
          str += String.fromCharCode(u0);
        } else {
          var ch = u0 - 0x10000;
          str +=
              String.fromCharCode(0xD800 | (ch >> 10), 0xDC00 | (ch & 0x3FF));
        }
      }
      return str;
    };

    /**
     * Given a pointer 'ptr' to a null-terminated UTF8-encoded string in the
     * emscripten HEAP, returns a copy of that string as a Javascript String
     * object.
     *
     * @param {number} ptr
     * @param {number=} maxBytesToRead - An optional length that specifies the
     *   maximum number of bytes to read. You can omit this parameter to scan
     * the string until the first 0 byte. If maxBytesToRead is passed, and the
     * string at [ptr, ptr+maxBytesToReadr[ contains a null byte in the middle,
     * then the string will cut short at that byte index.
     * @param {boolean=} ignoreNul - If true, the function will not stop on a
     *     NUL character.
     * @return {string}
     */
    var UTF8ToString = (ptr, maxBytesToRead, ignoreNul) => {
      assert(typeof ptr == 'number',
             `UTF8ToString expects a number (got ${typeof ptr})`);
      return ptr ? UTF8ArrayToString(HEAPU8, ptr, maxBytesToRead, ignoreNul)
                 : '';
    };
    var ___assert_fail = (condition, filename, line, func) =>
        abort(`Assertion failed: ${UTF8ToString(condition)}, at: ` + [
          filename ? UTF8ToString(filename) : 'unknown filename', line,
          func ? UTF8ToString(func) : 'unknown function'
        ]);

    var wasmTableMirror = [];

    /** @type {WebAssembly.Table} */
    var wasmTable;
    var getWasmTableEntry = (funcPtr) => {
      var func = wasmTableMirror[funcPtr];
      if (!func) {
        /** @suppress {checkTypes} */
        wasmTableMirror[funcPtr] = func = wasmTable.get(funcPtr);
      }
      /** @suppress {checkTypes} */
      assert(wasmTable.get(funcPtr) == func,
             'JavaScript-side Wasm function table mirror is out of date!');
      return func;
    };
    var ___call_sighandler = (fp, sig) => getWasmTableEntry(fp)(sig);

    /** @suppress {duplicate } */
    var syscallGetVarargI = () => {
      assert(SYSCALLS.varargs != undefined);
      // the `+` prepended here is necessary to convince the JSCompiler that
      // varargs is indeed a number.
      var ret = HEAP32[((+SYSCALLS.varargs) >> 2)];
      SYSCALLS.varargs += 4;
      return ret;
    };
    var syscallGetVarargP = syscallGetVarargI;

    var PATH = {
      isAbs : (path) => path.charAt(0) === '/',
      splitPath : (filename) => {
        var splitPathRe =
            /^(\/?|)([\s\S]*?)((?:\.{1,2}|[^\/]+?|)(\.[^.\/]*|))(?:[\/]*)$/;
        return splitPathRe.exec(filename).slice(1);
      },
      normalizeArray : (parts, allowAboveRoot) => {
        // if the path tries to go above the root, `up` ends up > 0
        var up = 0;
        for (var i = parts.length - 1; i >= 0; i--) {
          var last = parts[i];
          if (last === '.') {
            parts.splice(i, 1);
          } else if (last === '..') {
            parts.splice(i, 1);
            up++;
          } else if (up) {
            parts.splice(i, 1);
            up--;
          }
        }
        // if the path is allowed to go above the root, restore leading ..s
        if (allowAboveRoot) {
          for (; up; up--) {
            parts.unshift('..');
          }
        }
        return parts;
      },
      normalize : (path) => {
        var isAbsolute = PATH.isAbs(path),
            trailingSlash = path.slice(-1) === '/';
        // Normalize the path
        path =
            PATH.normalizeArray(path.split('/').filter((p) => !!p), !isAbsolute)
                .join('/');
        if (!path && !isAbsolute) {
          path = '.';
        }
        if (path && trailingSlash) {
          path += '/';
        }
        return (isAbsolute ? '/' : '') + path;
      },
      dirname : (path) => {
        var result = PATH.splitPath(path), root = result[0], dir = result[1];
        if (!root && !dir) {
          // No dirname whatsoever
          return '.';
        }
        if (dir) {
          // It has a dirname, strip trailing slash
          dir = dir.slice(0, -1);
        }
        return root + dir;
      },
      basename : (path) => path && path.match(/([^\/]+|\/)\/*$/)[1],
      join : (...paths) => PATH.normalize(paths.join('/')),
      join2 : (l, r) => PATH.normalize(l + '/' + r),
    };

    var initRandomFill =
        () => { return (view) => crypto.getRandomValues(view); };
    var randomFill = (view) => {
      // Lazily init on the first invocation.
      (randomFill = initRandomFill())(view);
    };

    var PATH_FS = {
      resolve : (...args) => {
        var resolvedPath = '', resolvedAbsolute = false;
        for (var i = args.length - 1; i >= -1 && !resolvedAbsolute; i--) {
          var path = (i >= 0) ? args[i] : FS.cwd();
          // Skip empty and invalid entries
          if (typeof path != 'string') {
            throw new TypeError('Arguments to path.resolve must be strings');
          } else if (!path) {
            return ''; // an invalid portion invalidates the whole thing
          }
          resolvedPath = path + '/' + resolvedPath;
          resolvedAbsolute = PATH.isAbs(path);
        }
        // At this point the path should be resolved to a full absolute path,
        // but handle relative paths to be safe (might happen when process.cwd()
        // fails)
        resolvedPath =
            PATH.normalizeArray(resolvedPath.split('/').filter((p) => !!p),
                                !resolvedAbsolute)
                .join('/');
        return ((resolvedAbsolute ? '/' : '') + resolvedPath) || '.';
      },
      relative : (from, to) => {
        from = PATH_FS.resolve(from).slice(1);
        to = PATH_FS.resolve(to).slice(1);
        function trim(arr) {
          var start = 0;
          for (; start < arr.length; start++) {
            if (arr[start] !== '')
              break;
          }
          var end = arr.length - 1;
          for (; end >= 0; end--) {
            if (arr[end] !== '')
              break;
          }
          if (start > end)
            return [];
          return arr.slice(start, end - start + 1);
        }
        var fromParts = trim(from.split('/'));
        var toParts = trim(to.split('/'));
        var length = Math.min(fromParts.length, toParts.length);
        var samePartsLength = length;
        for (var i = 0; i < length; i++) {
          if (fromParts[i] !== toParts[i]) {
            samePartsLength = i;
            break;
          }
        }
        var outputParts = [];
        for (var i = samePartsLength; i < fromParts.length; i++) {
          outputParts.push('..');
        }
        outputParts = outputParts.concat(toParts.slice(samePartsLength));
        return outputParts.join('/');
      },
    };

    var FS_stdin_getChar_buffer = [];

    var lengthBytesUTF8 = (str) => {
      var len = 0;
      for (var i = 0; i < str.length; ++i) {
        // Gotcha: charCodeAt returns a 16-bit word that is a UTF-16 encoded
        // code unit, not a Unicode code point of the character! So decode
        // UTF16->UTF32->UTF8.
        // See http://unicode.org/faq/utf_bom.html#utf16-3
        var c = str.charCodeAt(i); // possibly a lead surrogate
        if (c <= 0x7F) {
          len++;
        } else if (c <= 0x7FF) {
          len += 2;
        } else if (c >= 0xD800 && c <= 0xDFFF) {
          len += 4;
          ++i;
        } else {
          len += 3;
        }
      }
      return len;
    };

    var stringToUTF8Array = (str, heap, outIdx, maxBytesToWrite) => {
      assert(typeof str === 'string',
             `stringToUTF8Array expects a string (got ${typeof str})`);
      // Parameter maxBytesToWrite is not optional. Negative values, 0, null,
      // undefined and false each don't write out any bytes.
      if (!(maxBytesToWrite > 0))
        return 0;

      var startIdx = outIdx;
      var endIdx =
          outIdx + maxBytesToWrite - 1; // -1 for string null terminator.
      for (var i = 0; i < str.length; ++i) {
        // For UTF8 byte structure, see
        // http://en.wikipedia.org/wiki/UTF-8#Description and
        // https://www.ietf.org/rfc/rfc2279.txt and
        // https://tools.ietf.org/html/rfc3629
        var u = str.codePointAt(i);
        if (u <= 0x7F) {
          if (outIdx >= endIdx)
            break;
          heap[outIdx++] = u;
        } else if (u <= 0x7FF) {
          if (outIdx + 1 >= endIdx)
            break;
          heap[outIdx++] = 0xC0 | (u >> 6);
          heap[outIdx++] = 0x80 | (u & 63);
        } else if (u <= 0xFFFF) {
          if (outIdx + 2 >= endIdx)
            break;
          heap[outIdx++] = 0xE0 | (u >> 12);
          heap[outIdx++] = 0x80 | ((u >> 6) & 63);
          heap[outIdx++] = 0x80 | (u & 63);
        } else {
          if (outIdx + 3 >= endIdx)
            break;
          if (u > 0x10FFFF)
            warnOnce(
                'Invalid Unicode code point ' + ptrToString(u) +
                ' encountered when serializing a JS string to a UTF-8 string in wasm memory! (Valid unicode code points should be in range 0-0x10FFFF).');
          heap[outIdx++] = 0xF0 | (u >> 18);
          heap[outIdx++] = 0x80 | ((u >> 12) & 63);
          heap[outIdx++] = 0x80 | ((u >> 6) & 63);
          heap[outIdx++] = 0x80 | (u & 63);
          // Gotcha: if codePoint is over 0xFFFF, it is represented as a
          // surrogate pair in UTF-16. We need to manually skip over the second
          // code unit for correct iteration.
          i++;
        }
      }
      // Null-terminate the pointer to the buffer.
      heap[outIdx] = 0;
      return outIdx - startIdx;
    };
    /** @type {function(string, boolean=, number=)} */
    var intArrayFromString = (stringy, dontAddNull, length) => {
      var len = length > 0 ? length : lengthBytesUTF8(stringy) + 1;
      var u8array = new Array(len);
      var numBytesWritten =
          stringToUTF8Array(stringy, u8array, 0, u8array.length);
      if (dontAddNull)
        u8array.length = numBytesWritten;
      return u8array;
    };
    var FS_stdin_getChar = () => {
      if (!FS_stdin_getChar_buffer.length) {
        var result = null;
        if (typeof window != 'undefined' &&
            typeof window.prompt == 'function') {
          // Browser.
          result = window.prompt('Input: '); // returns null on cancel
          if (result !== null) {
            result += '\n';
          }
        } else {
        }
        if (!result) {
          return null;
        }
        FS_stdin_getChar_buffer = intArrayFromString(result, true);
      }
      return FS_stdin_getChar_buffer.shift();
    };
    var TTY = {
      ttys : [],
      init() {
        // https://github.com/emscripten-core/emscripten/pull/1555
        // if (ENVIRONMENT_IS_NODE) {
        //   // currently, FS.init does not distinguish if process.stdin is a
        //   file or TTY
        //   // device, it always assumes it's a TTY device. because of this,
        //   we're forcing
        //   // process.stdin to UTF8 encoding to at least make stdin reading
        //   compatible
        //   // with text files until FS.init can be refactored.
        //   process.stdin.setEncoding('utf8');
        // }
      },
      shutdown() {
        // https://github.com/emscripten-core/emscripten/pull/1555
        // if (ENVIRONMENT_IS_NODE) {
        //   // inolen: any idea as to why node -e 'process.stdin.read()'
        //   wouldn't exit immediately (with process.stdin being a tty)?
        //   // isaacs: because now it's reading from the stream, you've
        //   expressed interest in it, so that read() kicks off a _read() which
        //   creates a ReadReq operation
        //   // inolen: I thought read() in that case was a synchronous
        //   operation that just grabbed some amount of buffered data if it
        //   exists?
        //   // isaacs: it is. but it also triggers a _read() call, which calls
        //   readStart() on the handle
        //   // isaacs: do process.stdin.pause() and i'd think it'd probably
        //   close the pending call process.stdin.pause();
        // }
      },
      register(dev, ops) {
        TTY.ttys[dev] = {input : [], output : [], ops : ops};
        FS.registerDevice(dev, TTY.stream_ops);
      },
      stream_ops : {
        open(stream) {
          var tty = TTY.ttys[stream.node.rdev];
          if (!tty) {
            throw new FS.ErrnoError(43);
          }
          stream.tty = tty;
          stream.seekable = false;
        },
        close(stream) {
          // flush any pending line data
          stream.tty.ops.fsync(stream.tty);
        },
        fsync(stream) { stream.tty.ops.fsync(stream.tty); },
        read(stream, buffer, offset, length, pos /* ignored */) {
          if (!stream.tty || !stream.tty.ops.get_char) {
            throw new FS.ErrnoError(60);
          }
          var bytesRead = 0;
          for (var i = 0; i < length; i++) {
            var result;
            try {
              result = stream.tty.ops.get_char(stream.tty);
            } catch (e) {
              throw new FS.ErrnoError(29);
            }
            if (result === undefined && bytesRead === 0) {
              throw new FS.ErrnoError(6);
            }
            if (result === null || result === undefined)
              break;
            bytesRead++;
            buffer[offset + i] = result;
          }
          if (bytesRead) {
            stream.node.atime = Date.now();
          }
          return bytesRead;
        },
        write(stream, buffer, offset, length, pos) {
          if (!stream.tty || !stream.tty.ops.put_char) {
            throw new FS.ErrnoError(60);
          }
          try {
            for (var i = 0; i < length; i++) {
              stream.tty.ops.put_char(stream.tty, buffer[offset + i]);
            }
          } catch (e) {
            throw new FS.ErrnoError(29);
          }
          if (length) {
            stream.node.mtime = stream.node.ctime = Date.now();
          }
          return i;
        },
      },
      default_tty_ops : {
        get_char(tty) { return FS_stdin_getChar(); },
        put_char(tty, val) {
          if (val === null || val === 10) {
            out(UTF8ArrayToString(tty.output));
            tty.output = [];
          } else {
            if (val != 0)
              tty.output.push(
                  val); // val == 0 would cut text output off in the middle.
          }
        },
        fsync(tty) {
          if (tty.output?.length > 0) {
            out(UTF8ArrayToString(tty.output));
            tty.output = [];
          }
        },
        ioctl_tcgets(tty) {
          // typical setting
          return {
            c_iflag : 25856,
            c_oflag : 5,
            c_cflag : 191,
            c_lflag : 35387,
            c_cc : [
              0x03, 0x1c, 0x7f, 0x15, 0x04, 0x00, 0x01, 0x00, 0x11, 0x13, 0x1a,
              0x00, 0x12, 0x0f, 0x17, 0x16, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
              0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            ]
          };
        },
        ioctl_tcsets(tty, optional_actions, data) {
          // currently just ignore
          return 0;
        },
        ioctl_tiocgwinsz(tty) { return [ 24, 80 ]; },
      },
      default_tty1_ops : {
        put_char(tty, val) {
          if (val === null || val === 10) {
            err(UTF8ArrayToString(tty.output));
            tty.output = [];
          } else {
            if (val != 0)
              tty.output.push(val);
          }
        },
        fsync(tty) {
          if (tty.output?.length > 0) {
            err(UTF8ArrayToString(tty.output));
            tty.output = [];
          }
        },
      },
    };

    var mmapAlloc = (size) => {
      abort(
          'internal error: mmapAlloc called but `emscripten_builtin_memalign` native symbol not exported');
    };
    var MEMFS = {
      ops_table : null,
      mount(mount) { return MEMFS.createNode(null, '/', 16895, 0); },
      createNode(parent, name, mode, dev) {
        if (FS.isBlkdev(mode) || FS.isFIFO(mode)) {
          // no supported
          throw new FS.ErrnoError(63);
        }
        MEMFS.ops_table ||= {
          dir : {
            node : {
              getattr : MEMFS.node_ops.getattr,
              setattr : MEMFS.node_ops.setattr,
              lookup : MEMFS.node_ops.lookup,
              mknod : MEMFS.node_ops.mknod,
              rename : MEMFS.node_ops.rename,
              unlink : MEMFS.node_ops.unlink,
              rmdir : MEMFS.node_ops.rmdir,
              readdir : MEMFS.node_ops.readdir,
              symlink : MEMFS.node_ops.symlink
            },
            stream : {llseek : MEMFS.stream_ops.llseek}
          },
          file : {
            node : {
              getattr : MEMFS.node_ops.getattr,
              setattr : MEMFS.node_ops.setattr
            },
            stream : {
              llseek : MEMFS.stream_ops.llseek,
              read : MEMFS.stream_ops.read,
              write : MEMFS.stream_ops.write,
              mmap : MEMFS.stream_ops.mmap,
              msync : MEMFS.stream_ops.msync
            }
          },
          link : {
            node : {
              getattr : MEMFS.node_ops.getattr,
              setattr : MEMFS.node_ops.setattr,
              readlink : MEMFS.node_ops.readlink
            },
            stream : {}
          },
          chrdev : {
            node : {
              getattr : MEMFS.node_ops.getattr,
              setattr : MEMFS.node_ops.setattr
            },
            stream : FS.chrdev_stream_ops
          }
        };
        var node = FS.createNode(parent, name, mode, dev);
        if (FS.isDir(node.mode)) {
          node.node_ops = MEMFS.ops_table.dir.node;
          node.stream_ops = MEMFS.ops_table.dir.stream;
          node.contents = {};
        } else if (FS.isFile(node.mode)) {
          node.node_ops = MEMFS.ops_table.file.node;
          node.stream_ops = MEMFS.ops_table.file.stream;
          node.usedBytes =
              0; // The actual number of bytes used in the typed array, as
                 // opposed to contents.length which gives the whole capacity.
          // When the byte data of the file is populated, this will point to
          // either a typed array, or a normal JS array. Typed arrays are
          // preferred for performance, and used by default. However, typed
          // arrays are not resizable like normal JS arrays are, so there is a
          // small disk size penalty involved for appending file writes that
          // continuously grow a file similar to std::vector capacity vs used
          // -scheme.
          node.contents = null;
        } else if (FS.isLink(node.mode)) {
          node.node_ops = MEMFS.ops_table.link.node;
          node.stream_ops = MEMFS.ops_table.link.stream;
        } else if (FS.isChrdev(node.mode)) {
          node.node_ops = MEMFS.ops_table.chrdev.node;
          node.stream_ops = MEMFS.ops_table.chrdev.stream;
        }
        node.atime = node.mtime = node.ctime = Date.now();
        // add the new node to the parent
        if (parent) {
          parent.contents[name] = node;
          parent.atime = parent.mtime = parent.ctime = node.atime;
        }
        return node;
      },
      getFileDataAsTypedArray(node) {
        if (!node.contents)
          return new Uint8Array(0);
        if (node.contents.subarray)
          return node.contents.subarray(
              0,
              node.usedBytes); // Make sure to not return excess unused bytes.
        return new Uint8Array(node.contents);
      },
      expandFileStorage(node, newCapacity) {
        var prevCapacity = node.contents ? node.contents.length : 0;
        if (prevCapacity >= newCapacity)
          return; // No need to expand, the storage was already large enough.
        // Don't expand strictly to the given requested limit if it's only a
        // very small increase, but instead geometrically grow capacity. For
        // small filesizes (<1MB), perform size*2 geometric increase, but for
        // large sizes, do a much more conservative size*1.125 increase to avoid
        // overshooting the allocation cap by a very large margin.
        var CAPACITY_DOUBLING_MAX = 1024 * 1024;
        newCapacity =
            Math.max(newCapacity,
                     (prevCapacity *
                      (prevCapacity < CAPACITY_DOUBLING_MAX ? 2.0 : 1.125)) >>>
                         0);
        if (prevCapacity != 0)
          newCapacity = Math.max(
              newCapacity,
              256); // At minimum allocate 256b for each file when expanding.
        var oldContents = node.contents;
        node.contents = new Uint8Array(newCapacity); // Allocate new storage.
        if (node.usedBytes > 0)
          node.contents.set(oldContents.subarray(0, node.usedBytes),
                            0); // Copy old data over to the new storage.
      },
      resizeFileStorage(node, newSize) {
        if (node.usedBytes == newSize)
          return;
        if (newSize == 0) {
          node.contents =
              null; // Fully decommit when requesting a resize to zero.
          node.usedBytes = 0;
        } else {
          var oldContents = node.contents;
          node.contents = new Uint8Array(newSize); // Allocate new storage.
          if (oldContents) {
            node.contents.set(oldContents.subarray(
                0,
                Math.min(
                    newSize,
                    node.usedBytes))); // Copy old data over to the new storage.
          }
          node.usedBytes = newSize;
        }
      },
      node_ops : {
        getattr(node) {
          var attr = {};
          // device numbers reuse inode numbers.
          attr.dev = FS.isChrdev(node.mode) ? node.id : 1;
          attr.ino = node.id;
          attr.mode = node.mode;
          attr.nlink = 1;
          attr.uid = 0;
          attr.gid = 0;
          attr.rdev = node.rdev;
          if (FS.isDir(node.mode)) {
            attr.size = 4096;
          } else if (FS.isFile(node.mode)) {
            attr.size = node.usedBytes;
          } else if (FS.isLink(node.mode)) {
            attr.size = node.link.length;
          } else {
            attr.size = 0;
          }
          attr.atime = new Date(node.atime);
          attr.mtime = new Date(node.mtime);
          attr.ctime = new Date(node.ctime);
          // NOTE: In our implementation, st_blocks =
          // Math.ceil(st_size/st_blksize),
          //       but this is not required by the standard.
          attr.blksize = 4096;
          attr.blocks = Math.ceil(attr.size / attr.blksize);
          return attr;
        },
        setattr(node, attr) {
          for (const key of ["mode", "atime", "mtime", "ctime"]) {
            if (attr[key] != null) {
              node[key] = attr[key];
            }
          }
          if (attr.size !== undefined) {
            MEMFS.resizeFileStorage(node, attr.size);
          }
        },
        lookup(parent, name) { throw new FS.ErrnoError(44); },
        mknod(parent, name, mode,
              dev) { return MEMFS.createNode(parent, name, mode, dev); },
        rename(old_node, new_dir, new_name) {
          var new_node;
          try {
            new_node = FS.lookupNode(new_dir, new_name);
          } catch (e) {
          }
          if (new_node) {
            if (FS.isDir(old_node.mode)) {
              // if we're overwriting a directory at new_name, make sure it's
              // empty.
              for (var i in new_node.contents) {
                throw new FS.ErrnoError(55);
              }
            }
            FS.hashRemoveNode(new_node);
          }
          // do the internal rewiring
          delete old_node.parent.contents[old_node.name];
          new_dir.contents[new_name] = old_node;
          old_node.name = new_name;
          new_dir.ctime = new_dir.mtime = old_node.parent.ctime =
              old_node.parent.mtime = Date.now();
        },
        unlink(parent, name) {
          delete parent.contents[name];
          parent.ctime = parent.mtime = Date.now();
        },
        rmdir(parent, name) {
          var node = FS.lookupNode(parent, name);
          for (var i in node.contents) {
            throw new FS.ErrnoError(55);
          }
          delete parent.contents[name];
          parent.ctime = parent.mtime = Date.now();
        },
        readdir(node) { return [ '.', '..', ...Object.keys(node.contents) ]; },
        symlink(parent, newname, oldpath) {
          var node = MEMFS.createNode(parent, newname, 0o777 | 40960, 0);
          node.link = oldpath;
          return node;
        },
        readlink(node) {
          if (!FS.isLink(node.mode)) {
            throw new FS.ErrnoError(28);
          }
          return node.link;
        },
      },
      stream_ops : {
        read(stream, buffer, offset, length, position) {
          var contents = stream.node.contents;
          if (position >= stream.node.usedBytes)
            return 0;
          var size = Math.min(stream.node.usedBytes - position, length);
          assert(size >= 0);
          if (size > 8 && contents.subarray) { // non-trivial, and typed array
            buffer.set(contents.subarray(position, position + size), offset);
          } else {
            for (var i = 0; i < size; i++)
              buffer[offset + i] = contents[position + i];
          }
          return size;
        },
        write(stream, buffer, offset, length, position, canOwn) {
          // The data buffer should be a typed array view
          assert(!(buffer instanceof ArrayBuffer));
          // If the buffer is located in main memory (HEAP), and if
          // memory can grow, we can't hold on to references of the
          // memory buffer, as they may get invalidated. That means we
          // need to do copy its contents.
          if (buffer.buffer === HEAP8.buffer) {
            canOwn = false;
          }

          if (!length)
            return 0;
          var node = stream.node;
          node.mtime = node.ctime = Date.now();

          if (buffer.subarray &&
              (!node.contents ||
               node.contents.subarray)) { // This write is from a typed array to
                                          // a typed array?
            if (canOwn) {
              assert(position === 0,
                     'canOwn must imply no weird position inside the file');
              node.contents = buffer.subarray(offset, offset + length);
              node.usedBytes = length;
              return length;
            } else if (node.usedBytes === 0 &&
                       position ===
                           0) { // If this is a simple first write to an empty
                                // file, do a fast set since we don't need to
                                // care about old data.
              node.contents = buffer.slice(offset, offset + length);
              node.usedBytes = length;
              return length;
            } else if (position + length <=
                       node.usedBytes) { // Writing to an already allocated and
                                         // used subrange of the file?
              node.contents.set(buffer.subarray(offset, offset + length),
                                position);
              return length;
            }
          }

          // Appending to an existing file and we need to reallocate, or source
          // data did not come as a typed array.
          MEMFS.expandFileStorage(node, position + length);
          if (node.contents.subarray && buffer.subarray) {
            // Use typed array write which is available.
            node.contents.set(buffer.subarray(offset, offset + length),
                              position);
          } else {
            for (var i = 0; i < length; i++) {
              node.contents[position + i] =
                  buffer[offset + i]; // Or fall back to manual write if not.
            }
          }
          node.usedBytes = Math.max(node.usedBytes, position + length);
          return length;
        },
        llseek(stream, offset, whence) {
          var position = offset;
          if (whence === 1) {
            position += stream.position;
          } else if (whence === 2) {
            if (FS.isFile(stream.node.mode)) {
              position += stream.node.usedBytes;
            }
          }
          if (position < 0) {
            throw new FS.ErrnoError(28);
          }
          return position;
        },
        mmap(stream, length, position, prot, flags) {
          if (!FS.isFile(stream.node.mode)) {
            throw new FS.ErrnoError(43);
          }
          var ptr;
          var allocated;
          var contents = stream.node.contents;
          // Only make a new copy when MAP_PRIVATE is specified.
          if (!(flags & 2) && contents && contents.buffer === HEAP8.buffer) {
            // We can't emulate MAP_SHARED when the file is not backed by the
            // buffer we're mapping to (e.g. the HEAP buffer).
            allocated = false;
            ptr = contents.byteOffset;
          } else {
            allocated = true;
            ptr = mmapAlloc(length);
            if (!ptr) {
              throw new FS.ErrnoError(48);
            }
            if (contents) {
              // Try to avoid unnecessary slices.
              if (position > 0 || position + length < contents.length) {
                if (contents.subarray) {
                  contents = contents.subarray(position, position + length);
                } else {
                  contents = Array.prototype.slice.call(contents, position,
                                                        position + length);
                }
              }
              HEAP8.set(contents, ptr);
            }
          }
          return {ptr, allocated};
        },
        msync(stream, buffer, offset, length, mmapFlags) {
          MEMFS.stream_ops.write(stream, buffer, 0, length, offset, false);
          // should we check if bytesWritten and length are the same?
          return 0;
        },
      },
    };

    var FS_modeStringToFlags = (str) => {
      var flagModes = {
        'r' : 0,
        'r+' : 2,
        'w' : 512 | 64 | 1,
        'w+' : 512 | 64 | 2,
        'a' : 1024 | 64 | 1,
        'a+' : 1024 | 64 | 2,
      };
      var flags = flagModes[str];
      if (typeof flags == 'undefined') {
        throw new Error(`Unknown file open mode: ${str}`);
      }
      return flags;
    };

    var FS_getMode = (canRead, canWrite) => {
      var mode = 0;
      if (canRead)
        mode |= 292 | 73;
      if (canWrite)
        mode |= 146;
      return mode;
    };

    var WORKERFS = {
      DIR_MODE : 16895,
      FILE_MODE : 33279,
      reader : null,
      mount(mount) {
        assert(ENVIRONMENT_IS_WORKER);
        WORKERFS.reader ??= new FileReaderSync();
        var root = WORKERFS.createNode(null, '/', WORKERFS.DIR_MODE, 0);
        var createdParents = {};
        function ensureParent(path) {
          // return the parent node, creating subdirs as necessary
          var parts = path.split('/');
          var parent = root;
          for (var i = 0; i < parts.length - 1; i++) {
            var curr = parts.slice(0, i + 1).join('/');
            // Issue 4254: Using curr as a node name will prevent the node
            // from being found in FS.nameTable when FS.open is called on
            // a path which holds a child of this node,
            // given that all FS functions assume node names
            // are just their corresponding parts within their given path,
            // rather than incremental aggregates which include their parent's
            // directories.
            createdParents[curr] ||=
                WORKERFS.createNode(parent, parts[i], WORKERFS.DIR_MODE, 0);
            parent = createdParents[curr];
          }
          return parent;
        }
        function base(path) {
          var parts = path.split('/');
          return parts[parts.length - 1];
        }
        // We also accept FileList here, by using Array.prototype
        Array.prototype.forEach.call(mount.opts["files"] || [], function(file) {
          WORKERFS.createNode(ensureParent(file.name), base(file.name),
                              WORKERFS.FILE_MODE, 0, file,
                              file.lastModifiedDate);
        });
        (mount.opts["blobs"] || []).forEach((obj) => {
          WORKERFS.createNode(ensureParent(obj["name"]), base(obj["name"]),
                              WORKERFS.FILE_MODE, 0, obj["data"]);
        });
        (mount.opts["packages"] || []).forEach((pack) => {
          pack['metadata'].files.forEach((file) => {
            var name = file.filename.slice(1); // remove initial slash
            WORKERFS.createNode(ensureParent(name), base(name),
                                WORKERFS.FILE_MODE, 0,
                                pack['blob'].slice(file.start, file.end));
          });
        });
        return root;
      },
      createNode(parent, name, mode, dev, contents, mtime) {
        var node = FS.createNode(parent, name, mode);
        node.mode = mode;
        node.node_ops = WORKERFS.node_ops;
        node.stream_ops = WORKERFS.stream_ops;
        node.atime = node.mtime = node.ctime = (mtime || new Date).getTime();
        assert(WORKERFS.FILE_MODE !== WORKERFS.DIR_MODE);
        if (mode === WORKERFS.FILE_MODE) {
          node.size = contents.size;
          node.contents = contents;
        } else {
          node.size = 4096;
          node.contents = {};
        }
        if (parent) {
          parent.contents[name] = node;
        }
        return node;
      },
      node_ops : {
        getattr(node) {
          return {
            dev : 1,
            ino : node.id,
            mode : node.mode,
            nlink : 1,
            uid : 0,
            gid : 0,
            rdev : 0,
            size : node.size,
            atime : new Date(node.atime),
            mtime : new Date(node.mtime),
            ctime : new Date(node.ctime),
            blksize : 4096,
            blocks : Math.ceil(node.size / 4096),
          };
        },
        setattr(node, attr) {
          for (const key of ["mode", "atime", "mtime", "ctime"]) {
            if (attr[key] != null) {
              node[key] = attr[key];
            }
          }
        },
        lookup(parent, name) { throw new FS.ErrnoError(44); },
        mknod(parent, name, mode, dev) { throw new FS.ErrnoError(63); },
        rename(oldNode, newDir, newName) { throw new FS.ErrnoError(63); },
        unlink(parent, name) { throw new FS.ErrnoError(63); },
        rmdir(parent, name) { throw new FS.ErrnoError(63); },
        readdir(node) {
          var entries = [ '.', '..' ];
          for (var key of Object.keys(node.contents)) {
            entries.push(key);
          }
          return entries;
        },
        symlink(parent, newName, oldPath) { throw new FS.ErrnoError(63); },
      },
      stream_ops : {
        read(stream, buffer, offset, length, position) {
          if (position >= stream.node.size)
            return 0;
          var chunk = stream.node.contents.slice(position, position + length);
          var ab = WORKERFS.reader.readAsArrayBuffer(chunk);
          buffer.set(new Uint8Array(ab), offset);
          return chunk.size;
        },
        write(stream, buffer, offset, length,
              position) { throw new FS.ErrnoError(29); },
        llseek(stream, offset, whence) {
          var position = offset;
          if (whence === 1) {
            position += stream.position;
          } else if (whence === 2) {
            if (FS.isFile(stream.node.mode)) {
              position += stream.node.size;
            }
          }
          if (position < 0) {
            throw new FS.ErrnoError(28);
          }
          return position;
        },
      },
    };

    var ERRNO_CODES = {
      'EPERM' : 63,
      'ENOENT' : 44,
      'ESRCH' : 71,
      'EINTR' : 27,
      'EIO' : 29,
      'ENXIO' : 60,
      'E2BIG' : 1,
      'ENOEXEC' : 45,
      'EBADF' : 8,
      'ECHILD' : 12,
      'EAGAIN' : 6,
      'EWOULDBLOCK' : 6,
      'ENOMEM' : 48,
      'EACCES' : 2,
      'EFAULT' : 21,
      'ENOTBLK' : 105,
      'EBUSY' : 10,
      'EEXIST' : 20,
      'EXDEV' : 75,
      'ENODEV' : 43,
      'ENOTDIR' : 54,
      'EISDIR' : 31,
      'EINVAL' : 28,
      'ENFILE' : 41,
      'EMFILE' : 33,
      'ENOTTY' : 59,
      'ETXTBSY' : 74,
      'EFBIG' : 22,
      'ENOSPC' : 51,
      'ESPIPE' : 70,
      'EROFS' : 69,
      'EMLINK' : 34,
      'EPIPE' : 64,
      'EDOM' : 18,
      'ERANGE' : 68,
      'ENOMSG' : 49,
      'EIDRM' : 24,
      'ECHRNG' : 106,
      'EL2NSYNC' : 156,
      'EL3HLT' : 107,
      'EL3RST' : 108,
      'ELNRNG' : 109,
      'EUNATCH' : 110,
      'ENOCSI' : 111,
      'EL2HLT' : 112,
      'EDEADLK' : 16,
      'ENOLCK' : 46,
      'EBADE' : 113,
      'EBADR' : 114,
      'EXFULL' : 115,
      'ENOANO' : 104,
      'EBADRQC' : 103,
      'EBADSLT' : 102,
      'EDEADLOCK' : 16,
      'EBFONT' : 101,
      'ENOSTR' : 100,
      'ENODATA' : 116,
      'ETIME' : 117,
      'ENOSR' : 118,
      'ENONET' : 119,
      'ENOPKG' : 120,
      'EREMOTE' : 121,
      'ENOLINK' : 47,
      'EADV' : 122,
      'ESRMNT' : 123,
      'ECOMM' : 124,
      'EPROTO' : 65,
      'EMULTIHOP' : 36,
      'EDOTDOT' : 125,
      'EBADMSG' : 9,
      'ENOTUNIQ' : 126,
      'EBADFD' : 127,
      'EREMCHG' : 128,
      'ELIBACC' : 129,
      'ELIBBAD' : 130,
      'ELIBSCN' : 131,
      'ELIBMAX' : 132,
      'ELIBEXEC' : 133,
      'ENOSYS' : 52,
      'ENOTEMPTY' : 55,
      'ENAMETOOLONG' : 37,
      'ELOOP' : 32,
      'EOPNOTSUPP' : 138,
      'EPFNOSUPPORT' : 139,
      'ECONNRESET' : 15,
      'ENOBUFS' : 42,
      'EAFNOSUPPORT' : 5,
      'EPROTOTYPE' : 67,
      'ENOTSOCK' : 57,
      'ENOPROTOOPT' : 50,
      'ESHUTDOWN' : 140,
      'ECONNREFUSED' : 14,
      'EADDRINUSE' : 3,
      'ECONNABORTED' : 13,
      'ENETUNREACH' : 40,
      'ENETDOWN' : 38,
      'ETIMEDOUT' : 73,
      'EHOSTDOWN' : 142,
      'EHOSTUNREACH' : 23,
      'EINPROGRESS' : 26,
      'EALREADY' : 7,
      'EDESTADDRREQ' : 17,
      'EMSGSIZE' : 35,
      'EPROTONOSUPPORT' : 66,
      'ESOCKTNOSUPPORT' : 137,
      'EADDRNOTAVAIL' : 4,
      'ENETRESET' : 39,
      'EISCONN' : 30,
      'ENOTCONN' : 53,
      'ETOOMANYREFS' : 141,
      'EUSERS' : 136,
      'EDQUOT' : 19,
      'ESTALE' : 72,
      'ENOTSUP' : 138,
      'ENOMEDIUM' : 148,
      'EILSEQ' : 25,
      'EOVERFLOW' : 61,
      'ECANCELED' : 11,
      'ENOTRECOVERABLE' : 56,
      'EOWNERDEAD' : 62,
      'ESTRPIPE' : 135,
    };
    var PROXYFS = {
      mount(mount) {
        return PROXYFS.createNode(null, '/',
                                  mount.opts.fs.lstat(mount.opts.root).mode, 0);
      },
      createNode(parent, name, mode, dev) {
        if (!FS.isDir(mode) && !FS.isFile(mode) && !FS.isLink(mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var node = FS.createNode(parent, name, mode);
        node.node_ops = PROXYFS.node_ops;
        node.stream_ops = PROXYFS.stream_ops;
        return node;
      },
      realPath(node) {
        var parts = [];
        while (node.parent !== node) {
          parts.push(node.name);
          node = node.parent;
        }
        parts.push(node.mount.opts.root);
        parts.reverse();
        return PATH.join(...parts);
      },
      node_ops : {
        getattr(node) {
          var path = PROXYFS.realPath(node);
          var stat;
          try {
            stat = node.mount.opts.fs.lstat(path);
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
          return {
            dev : stat.dev,
            ino : stat.ino,
            mode : stat.mode,
            nlink : stat.nlink,
            uid : stat.uid,
            gid : stat.gid,
            rdev : stat.rdev,
            size : stat.size,
            atime : stat.atime,
            mtime : stat.mtime,
            ctime : stat.ctime,
            blksize : stat.blksize,
            blocks : stat.blocks
          };
        },
        setattr(node, attr) {
          var path = PROXYFS.realPath(node);
          try {
            if (attr.mode !== undefined) {
              node.mount.opts.fs.chmod(path, attr.mode);
              // update the common node structure mode as well
              node.mode = attr.mode;
            }
            if (attr.atime || attr.mtime) {
              var atime = new Date(attr.atime || attr.mtime);
              var mtime = new Date(attr.mtime || attr.atime);
              node.mount.opts.fs.utime(path, atime, mtime);
            }
            if (attr.size !== undefined) {
              node.mount.opts.fs.truncate(path, attr.size);
            }
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        lookup(parent, name) {
          try {
            var path = PATH.join2(PROXYFS.realPath(parent), name);
            var mode = parent.mount.opts.fs.lstat(path).mode;
            var node = PROXYFS.createNode(parent, name, mode);
            return node;
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        mknod(parent, name, mode, dev) {
          var node = PROXYFS.createNode(parent, name, mode, dev);
          // create the backing node for this in the fs root as well
          var path = PROXYFS.realPath(node);
          try {
            if (FS.isDir(node.mode)) {
              node.mount.opts.fs.mkdir(path, node.mode);
            } else {
              node.mount.opts.fs.writeFile(path, '', {mode : node.mode});
            }
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
          return node;
        },
        rename(oldNode, newDir, newName) {
          var oldPath = PROXYFS.realPath(oldNode);
          var newPath = PATH.join2(PROXYFS.realPath(newDir), newName);
          try {
            oldNode.mount.opts.fs.rename(oldPath, newPath);
            oldNode.name = newName;
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        unlink(parent, name) {
          var path = PATH.join2(PROXYFS.realPath(parent), name);
          try {
            parent.mount.opts.fs.unlink(path);
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        rmdir(parent, name) {
          var path = PATH.join2(PROXYFS.realPath(parent), name);
          try {
            parent.mount.opts.fs.rmdir(path);
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        readdir(node) {
          var path = PROXYFS.realPath(node);
          try {
            return node.mount.opts.fs.readdir(path);
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        symlink(parent, newName, oldPath) {
          var newPath = PATH.join2(PROXYFS.realPath(parent), newName);
          try {
            parent.mount.opts.fs.symlink(oldPath, newPath);
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        readlink(node) {
          var path = PROXYFS.realPath(node);
          try {
            return node.mount.opts.fs.readlink(path);
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
      },
      stream_ops : {
        open(stream) {
          var path = PROXYFS.realPath(stream.node);
          try {
            stream.nfd = stream.node.mount.opts.fs.open(path, stream.flags);
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        close(stream) {
          try {
            stream.node.mount.opts.fs.close(stream.nfd);
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        read(stream, buffer, offset, length, position) {
          try {
            return stream.node.mount.opts.fs.read(stream.nfd, buffer, offset,
                                                  length, position);
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        write(stream, buffer, offset, length, position) {
          try {
            return stream.node.mount.opts.fs.write(stream.nfd, buffer, offset,
                                                   length, position);
          } catch (e) {
            if (!e.code)
              throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },
        llseek(stream, offset, whence) {
          var position = offset;
          if (whence === 1) {
            position += stream.position;
          } else if (whence === 2) {
            if (FS.isFile(stream.node.mode)) {
              try {
                var stat = stream.node.node_ops.getattr(stream.node);
                position += stat.size;
              } catch (e) {
                throw new FS.ErrnoError(ERRNO_CODES[e.code]);
              }
            }
          }

          if (position < 0) {
            throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
          }

          return position;
        },
      },
    };

    var strError = (errno) => UTF8ToString(_strerror(errno));

    var asyncLoad = async (url) => {
      var arrayBuffer = await readAsync(url);
      assert(arrayBuffer,
             `Loading data file "${url}" failed (no arrayBuffer).`);
      return new Uint8Array(arrayBuffer);
    };

    var FS_createDataFile = (...args) => FS.createDataFile(...args);

    var getUniqueRunDependency = (id) => {
      var orig = id;
      while (1) {
        if (!runDependencyTracking[id])
          return id;
        id = orig + Math.random();
      }
    };

    var preloadPlugins = [];
    var FS_handledByPreloadPlugin = async (byteArray, fullname) => {
      // Ensure plugins are ready.
      if (typeof Browser != 'undefined')
        Browser.init();

      for (var plugin of preloadPlugins) {
        if (plugin['canHandle'](fullname)) {
          assert(
              plugin['handle'].constructor.name === 'AsyncFunction',
              'Filesystem plugin handlers must be async functions (See #24914)')
          return plugin['handle'](byteArray, fullname);
        }
      }
      // In no plugin handled this file then return the original/unmodified
      // byteArray.
      return byteArray;
    };
    var FS_preloadFile = async (parent, name, url, canRead, canWrite,
                                dontCreateFile, canOwn, preFinish) => {
      // TODO we should allow people to just pass in a complete filename instead
      // of parent and name being that we just join them anyways
      var fullname = name ? PATH_FS.resolve(PATH.join2(parent, name)) : parent;
      var dep = getUniqueRunDependency(
          `cp ${fullname}`); // might have several active requests for the same
                             // fullname
      addRunDependency(dep);

      try {
        var byteArray = url;
        if (typeof url == 'string') {
          byteArray = await asyncLoad(url);
        }

        byteArray = await FS_handledByPreloadPlugin(byteArray, fullname);
        preFinish?.();
        if (!dontCreateFile) {
          FS_createDataFile(parent, name, byteArray, canRead, canWrite, canOwn);
        }
      } finally {
        removeRunDependency(dep);
      }
    };
    var FS_createPreloadedFile =
        (parent, name, url, canRead, canWrite, onload, onerror, dontCreateFile,
         canOwn, preFinish) => {
          FS_preloadFile(parent, name, url, canRead, canWrite, dontCreateFile,
                         canOwn, preFinish)
              .then(onload)
              .catch(onerror);
        };
    var FS = {
      root : null,
      mounts : [],
      devices : {},
      streams : [],
      nextInode : 1,
      nameTable : null,
      currentPath : "/",
      initialized : false,
      ignorePermissions : true,
      filesystems : null,
      syncFSRequests : 0,
      readFiles : {},
      ErrnoError : class extends
          Error {
            name = 'ErrnoError';
            // We set the `name` property to be able to identify `FS.ErrnoError`
            // - the `name` is a standard ECMA-262 property of error objects.
            // Kind of good to have it anyway.
            // - when using PROXYFS, an error can come from an underlying FS
            // as different FS objects have their own FS.ErrnoError each,
            // the test `err instanceof FS.ErrnoError` won't detect an error
            // coming from another filesystem, causing bugs. we'll use the
            // reliable test `err.name == "ErrnoError"` instead
            constructor(errno) {
              super(runtimeInitialized ? strError(errno) : '');
              this.errno = errno;
              for (var key in ERRNO_CODES) {
                if (ERRNO_CODES[key] === errno) {
                  this.code = key;
                  break;
                }
              }
            }
          },
          FSStream : class {
            shared = {};
            get object() { return this.node; }
            set object(val) { this.node = val; }
            get isRead() { return (this.flags & 2097155) !== 1; }
            get isWrite() { return (this.flags & 2097155) !== 0; }
            get isAppend() { return (this.flags & 1024); }
            get flags() { return this.shared.flags; }
            set flags(val) { this.shared.flags = val; }
            get position() { return this.shared.position; }
            set position(val) { this.shared.position = val; }
          },
          FSNode : class {
            node_ops = {};
            stream_ops = {};
            readMode = 292 | 73;
            writeMode = 146;
            mounted = null;
            constructor(parent, name, mode, rdev) {
              if (!parent) {
                parent = this; // root node sets parent to itself
              }
              this.parent = parent;
              this.mount = parent.mount;
              this.id = FS.nextInode++;
              this.name = name;
              this.mode = mode;
              this.rdev = rdev;
              this.atime = this.mtime = this.ctime = Date.now();
            }
            get read() { return (this.mode & this.readMode) === this.readMode; }
            set read(val) {
              val ? this.mode |= this.readMode : this.mode &= ~this.readMode;
            }
            get write() {
              return (this.mode & this.writeMode) === this.writeMode;
            }
            set write(val) {
              val ? this.mode |= this.writeMode : this.mode &= ~this.writeMode;
            }
            get isFolder() { return FS.isDir(this.mode); }
            get isDevice() { return FS.isChrdev(this.mode); }
          },
          lookupPath(path, opts = {}) {
            if (!path) {
              throw new FS.ErrnoError(44);
            }
            opts.follow_mount ??= true

            if (!PATH.isAbs(path)) {
              path = FS.cwd() + '/' + path;
            }

            // limit max consecutive symlinks to 40 (SYMLOOP_MAX).
            linkloop: for (var nlinks = 0; nlinks < 40; nlinks++) {
              // split the absolute path
              var parts = path.split('/').filter((p) => !!p);

              // start at the root
              var current = FS.root;
              var current_path = '/';

              for (var i = 0; i < parts.length; i++) {
                var islast = (i === parts.length - 1);
                if (islast && opts.parent) {
                  // stop resolving
                  break;
                }

                if (parts[i] === '.') {
                  continue;
                }

                if (parts[i] === '..') {
                  current_path = PATH.dirname(current_path);
                  if (FS.isRoot(current)) {
                    path = current_path + '/' + parts.slice(i + 1).join('/');
                    // We're making progress here, don't let many consecutive
                    // ..'s lead to ELOOP
                    nlinks--;
                    continue linkloop;
                  } else {
                    current = current.parent;
                  }
                  continue;
                }

                current_path = PATH.join2(current_path, parts[i]);
                try {
                  current = FS.lookupNode(current, parts[i]);
                } catch (e) {
                  // if noent_okay is true, suppress a ENOENT in the last
                  // component and return an object with an undefined node. This
                  // is needed for resolving symlinks in the path when creating
                  // a file.
                  if ((e?.errno === 44) && islast && opts.noent_okay) {
                    return {path : current_path};
                  }
                  throw e;
                }

                // jump to the mount's root node if this is a mountpoint
                if (FS.isMountpoint(current) &&
                    (!islast || opts.follow_mount)) {
                  current = current.mounted.root;
                }

                // by default, lookupPath will not follow a symlink if it is the
                // final path component. setting opts.follow = true will
                // override this behavior.
                if (FS.isLink(current.mode) && (!islast || opts.follow)) {
                  if (!current.node_ops.readlink) {
                    throw new FS.ErrnoError(52);
                  }
                  var link = current.node_ops.readlink(current);
                  if (!PATH.isAbs(link)) {
                    link = PATH.dirname(current_path) + '/' + link;
                  }
                  path = link + '/' + parts.slice(i + 1).join('/');
                  continue linkloop;
                }
              }
              return {path : current_path, node : current};
            }
            throw new FS.ErrnoError(32);
          },
          getPath(node) {
            var path;
            while (true) {
              if (FS.isRoot(node)) {
                var mount = node.mount.mountpoint;
                if (!path)
                  return mount;
                return mount[mount.length - 1] !== '/' ? `${mount}/${path}`
                                                       : mount + path;
              }
              path = path ? `${node.name}/${path}` : node.name;
              node = node.parent;
            }
          },
          hashName(parentid, name) {
            var hash = 0;

            for (var i = 0; i < name.length; i++) {
              hash = ((hash << 5) - hash + name.charCodeAt(i)) | 0;
            }
            return ((parentid + hash) >>> 0) % FS.nameTable.length;
          },
          hashAddNode(node) {
            var hash = FS.hashName(node.parent.id, node.name);
            node.name_next = FS.nameTable[hash];
            FS.nameTable[hash] = node;
          },
          hashRemoveNode(node) {
            var hash = FS.hashName(node.parent.id, node.name);
            if (FS.nameTable[hash] === node) {
              FS.nameTable[hash] = node.name_next;
            } else {
              var current = FS.nameTable[hash];
              while (current) {
                if (current.name_next === node) {
                  current.name_next = node.name_next;
                  break;
                }
                current = current.name_next;
              }
            }
          },
          lookupNode(parent, name) {
            var errCode = FS.mayLookup(parent);
            if (errCode) {
              throw new FS.ErrnoError(errCode);
            }
            var hash = FS.hashName(parent.id, name);
            for (var node = FS.nameTable[hash]; node; node = node.name_next) {
              var nodeName = node.name;
              if (node.parent.id === parent.id && nodeName === name) {
                return node;
              }
            }
            // if we failed to find it in the cache, call into the VFS
            return FS.lookup(parent, name);
          },
          createNode(parent, name, mode, rdev) {
            assert(typeof parent == 'object')
            var node = new FS.FSNode(parent, name, mode, rdev);

            FS.hashAddNode(node);

            return node;
          },
          destroyNode(node) { FS.hashRemoveNode(node); },
          isRoot(node) { return node === node.parent; },
          isMountpoint(node) { return !!node.mounted; },
          isFile(mode) { return (mode & 61440) === 32768; },
          isDir(mode) { return (mode & 61440) === 16384; },
          isLink(mode) { return (mode & 61440) === 40960; },
          isChrdev(mode) { return (mode & 61440) === 8192; },
          isBlkdev(mode) { return (mode & 61440) === 24576; },
          isFIFO(mode) { return (mode & 61440) === 4096; },
          isSocket(mode) { return (mode & 49152) === 49152; },
          flagsToPermissionString(flag) {
            var perms = [ 'r', 'w', 'rw' ][flag & 3];
            if ((flag & 512)) {
              perms += 'w';
            }
            return perms;
          },
          nodePermissions(node, perms) {
            if (FS.ignorePermissions) {
              return 0;
            }
            // return 0 if any user, group or owner bits are set.
            if (perms.includes('r') && !(node.mode & 292)) {
              return 2;
            } else if (perms.includes('w') && !(node.mode & 146)) {
              return 2;
            } else if (perms.includes('x') && !(node.mode & 73)) {
              return 2;
            }
            return 0;
          },
          mayLookup(dir) {
            if (!FS.isDir(dir.mode))
              return 54;
            var errCode = FS.nodePermissions(dir, 'x');
            if (errCode)
              return errCode;
            if (!dir.node_ops.lookup)
              return 2;
            return 0;
          },
          mayCreate(dir, name) {
            if (!FS.isDir(dir.mode)) {
              return 54;
            }
            try {
              var node = FS.lookupNode(dir, name);
              return 20;
            } catch (e) {
            }
            return FS.nodePermissions(dir, 'wx');
          },
          mayDelete(dir, name, isdir) {
            var node;
            try {
              node = FS.lookupNode(dir, name);
            } catch (e) {
              return e.errno;
            }
            var errCode = FS.nodePermissions(dir, 'wx');
            if (errCode) {
              return errCode;
            }
            if (isdir) {
              if (!FS.isDir(node.mode)) {
                return 54;
              }
              if (FS.isRoot(node) || FS.getPath(node) === FS.cwd()) {
                return 10;
              }
            } else {
              if (FS.isDir(node.mode)) {
                return 31;
              }
            }
            return 0;
          },
          mayOpen(node, flags) {
            if (!node) {
              return 44;
            }
            if (FS.isLink(node.mode)) {
              return 32;
            } else if (FS.isDir(node.mode)) {
              if (FS.flagsToPermissionString(flags) !== 'r' // opening for write
                  || (flags & (512 | 64))) { // TODO: check for O_SEARCH? (==
                                             // search for dir only)
                return 31;
              }
            }
            return FS.nodePermissions(node, FS.flagsToPermissionString(flags));
          },
          checkOpExists(op, err) {
            if (!op) {
              throw new FS.ErrnoError(err);
            }
            return op;
          },
          MAX_OPEN_FDS : 4096,
          nextfd() {
            for (var fd = 0; fd <= FS.MAX_OPEN_FDS; fd++) {
              if (!FS.streams[fd]) {
                return fd;
              }
            }
            throw new FS.ErrnoError(33);
          },
          getStreamChecked(fd) {
            var stream = FS.getStream(fd);
            if (!stream) {
              throw new FS.ErrnoError(8);
            }
            return stream;
          },
          getStream : (fd) => FS.streams[fd],
          createStream(stream, fd = -1) {
            assert(fd >= -1);

            // clone it, so we can return an instance of FSStream
            stream = Object.assign(new FS.FSStream(), stream);
            if (fd == -1) {
              fd = FS.nextfd();
            }
            stream.fd = fd;
            FS.streams[fd] = stream;
            return stream;
          },
          closeStream(fd) { FS.streams[fd] = null; },
          dupStream(origStream, fd = -1) {
            var stream = FS.createStream(origStream, fd);
            stream.stream_ops?.dup?.(stream);
            return stream;
          },
          doSetAttr(stream, node, attr) {
            var setattr = stream?.stream_ops.setattr;
            var arg = setattr ? stream : node;
            setattr ??= node.node_ops.setattr;
            FS.checkOpExists(setattr, 63)
            setattr(arg, attr);
          },
          chrdev_stream_ops : {
            open(stream) {
              var device = FS.getDevice(stream.node.rdev);
              // override node's stream ops with the device's
              stream.stream_ops = device.stream_ops;
              // forward the open call
              stream.stream_ops.open?.(stream);
            },
            llseek() { throw new FS.ErrnoError(70); },
          },
          major : (dev) => ((dev) >> 8),
          minor : (dev) => ((dev) & 0xff),
          makedev : (ma, mi) => ((ma) << 8 | (mi)),
          registerDevice(dev, ops) { FS.devices[dev] = {stream_ops : ops}; },
          getDevice : (dev) => FS.devices[dev],
          getMounts(mount) {
            var mounts = [];
            var check = [ mount ];

            while (check.length) {
              var m = check.pop();

              mounts.push(m);

              check.push(...m.mounts);
            }

            return mounts;
          },
          syncfs(populate, callback) {
            if (typeof populate == 'function') {
              callback = populate;
              populate = false;
            }

            FS.syncFSRequests++;

            if (FS.syncFSRequests > 1) {
              err(`warning: ${
                  FS.syncFSRequests} FS.syncfs operations in flight at once, probably just doing extra work`);
            }

            var mounts = FS.getMounts(FS.root.mount);
            var completed = 0;

            function doCallback(errCode) {
              assert(FS.syncFSRequests > 0);
              FS.syncFSRequests--;
              return callback(errCode);
            }

            function done(errCode) {
              if (errCode) {
                if (!done.errored) {
                  done.errored = true;
                  return doCallback(errCode);
                }
                return;
              }
              if (++completed >= mounts.length) {
                doCallback(null);
              }
            };

            // sync all mounts
            mounts.forEach((mount) => {
              if (!mount.type.syncfs) {
                return done(null);
              }
              mount.type.syncfs(mount, populate, done);
            });
          },
          mount(type, opts, mountpoint) {
            if (typeof type == 'string') {
              // The filesystem was not included, and instead we have an error
              // message stored in the variable.
              throw type;
            }
            var root = mountpoint === '/';
            var pseudo = !mountpoint;
            var node;

            if (root && FS.root) {
              throw new FS.ErrnoError(10);
            } else if (!root && !pseudo) {
              var lookup = FS.lookupPath(mountpoint, {follow_mount : false});

              mountpoint = lookup.path; // use the absolute path
              node = lookup.node;

              if (FS.isMountpoint(node)) {
                throw new FS.ErrnoError(10);
              }

              if (!FS.isDir(node.mode)) {
                throw new FS.ErrnoError(54);
              }
            }

            var mount = {type, opts, mountpoint, mounts : []};

            // create a root node for the fs
            var mountRoot = type.mount(mount);
            mountRoot.mount = mount;
            mount.root = mountRoot;

            if (root) {
              FS.root = mountRoot;
            } else if (node) {
              // set as a mountpoint
              node.mounted = mount;

              // add the new mount to the current mount's children
              if (node.mount) {
                node.mount.mounts.push(mount);
              }
            }

            return mountRoot;
          },
          unmount(mountpoint) {
            var lookup = FS.lookupPath(mountpoint, {follow_mount : false});

            if (!FS.isMountpoint(lookup.node)) {
              throw new FS.ErrnoError(28);
            }

            // destroy the nodes for this mount, and all its child mounts
            var node = lookup.node;
            var mount = node.mounted;
            var mounts = FS.getMounts(mount);

            Object.keys(FS.nameTable).forEach((hash) => {
              var current = FS.nameTable[hash];

              while (current) {
                var next = current.name_next;

                if (mounts.includes(current.mount)) {
                  FS.destroyNode(current);
                }

                current = next;
              }
            });

            // no longer a mountpoint
            node.mounted = null;

            // remove this mount from the child mounts
            var idx = node.mount.mounts.indexOf(mount);
            assert(idx !== -1);
            node.mount.mounts.splice(idx, 1);
          },
          lookup(parent, name) { return parent.node_ops.lookup(parent, name); },
          mknod(path, mode, dev) {
            var lookup = FS.lookupPath(path, {parent : true});
            var parent = lookup.node;
            var name = PATH.basename(path);
            if (!name) {
              throw new FS.ErrnoError(28);
            }
            if (name === '.' || name === '..') {
              throw new FS.ErrnoError(20);
            }
            var errCode = FS.mayCreate(parent, name);
            if (errCode) {
              throw new FS.ErrnoError(errCode);
            }
            if (!parent.node_ops.mknod) {
              throw new FS.ErrnoError(63);
            }
            return parent.node_ops.mknod(parent, name, mode, dev);
          },
          statfs(path) {
            return FS.statfsNode(FS.lookupPath(path, {follow : true}).node);
          },
          statfsStream(stream) {
            // We keep a separate statfsStream function because noderawfs
            // overrides it. In noderawfs, stream.node is sometimes null.
            // Instead, we need to look at stream.path.
            return FS.statfsNode(stream.node);
          },
          statfsNode(node) {
            // NOTE: None of the defaults here are true. We're just returning
            // safe and
            //       sane values. Currently nodefs and rawfs replace these
            //       defaults, other file systems leave them alone.
            var rtn = {
              bsize : 4096,
              frsize : 4096,
              blocks : 1e6,
              bfree : 5e5,
              bavail : 5e5,
              files : FS.nextInode,
              ffree : FS.nextInode - 1,
              fsid : 42,
              flags : 2,
              namelen : 255,
            };

            if (node.node_ops.statfs) {
              Object.assign(rtn, node.node_ops.statfs(node.mount.opts.root));
            }
            return rtn;
          },
          create(path, mode = 0o666) {
            mode &= 4095;
            mode |= 32768;
            return FS.mknod(path, mode, 0);
          },
          mkdir(path, mode = 0o777) {
            mode &= 511 | 512;
            mode |= 16384;
            return FS.mknod(path, mode, 0);
          },
          mkdirTree(path, mode) {
            var dirs = path.split('/');
            var d = '';
            for (var dir of dirs) {
              if (!dir)
                continue;
              if (d || PATH.isAbs(path))
                d += '/';
              d += dir;
              try {
                FS.mkdir(d, mode);
              } catch (e) {
                if (e.errno != 20)
                  throw e;
              }
            }
          },
          mkdev(path, mode, dev) {
            if (typeof dev == 'undefined') {
              dev = mode;
              mode = 0o666;
            }
            mode |= 8192;
            return FS.mknod(path, mode, dev);
          },
          symlink(oldpath, newpath) {
            if (!PATH_FS.resolve(oldpath)) {
              throw new FS.ErrnoError(44);
            }
            var lookup = FS.lookupPath(newpath, {parent : true});
            var parent = lookup.node;
            if (!parent) {
              throw new FS.ErrnoError(44);
            }
            var newname = PATH.basename(newpath);
            var errCode = FS.mayCreate(parent, newname);
            if (errCode) {
              throw new FS.ErrnoError(errCode);
            }
            if (!parent.node_ops.symlink) {
              throw new FS.ErrnoError(63);
            }
            return parent.node_ops.symlink(parent, newname, oldpath);
          },
          rename(old_path, new_path) {
            var old_dirname = PATH.dirname(old_path);
            var new_dirname = PATH.dirname(new_path);
            var old_name = PATH.basename(old_path);
            var new_name = PATH.basename(new_path);
            // parents must exist
            var lookup, old_dir, new_dir;

            // let the errors from non existent directories percolate up
            lookup = FS.lookupPath(old_path, {parent : true});
            old_dir = lookup.node;
            lookup = FS.lookupPath(new_path, {parent : true});
            new_dir = lookup.node;

            if (!old_dir || !new_dir)
              throw new FS.ErrnoError(44);
            // need to be part of the same mount
            if (old_dir.mount !== new_dir.mount) {
              throw new FS.ErrnoError(75);
            }
            // source must exist
            var old_node = FS.lookupNode(old_dir, old_name);
            // old path should not be an ancestor of the new path
            var relative = PATH_FS.relative(old_path, new_dirname);
            if (relative.charAt(0) !== '.') {
              throw new FS.ErrnoError(28);
            }
            // new path should not be an ancestor of the old path
            relative = PATH_FS.relative(new_path, old_dirname);
            if (relative.charAt(0) !== '.') {
              throw new FS.ErrnoError(55);
            }
            // see if the new path already exists
            var new_node;
            try {
              new_node = FS.lookupNode(new_dir, new_name);
            } catch (e) {
              // not fatal
            }
            // early out if nothing needs to change
            if (old_node === new_node) {
              return;
            }
            // we'll need to delete the old entry
            var isdir = FS.isDir(old_node.mode);
            var errCode = FS.mayDelete(old_dir, old_name, isdir);
            if (errCode) {
              throw new FS.ErrnoError(errCode);
            }
            // need delete permissions if we'll be overwriting.
            // need create permissions if new doesn't already exist.
            errCode = new_node ? FS.mayDelete(new_dir, new_name, isdir)
                               : FS.mayCreate(new_dir, new_name);
            if (errCode) {
              throw new FS.ErrnoError(errCode);
            }
            if (!old_dir.node_ops.rename) {
              throw new FS.ErrnoError(63);
            }
            if (FS.isMountpoint(old_node) ||
                (new_node && FS.isMountpoint(new_node))) {
              throw new FS.ErrnoError(10);
            }
            // if we are going to change the parent, check write permissions
            if (new_dir !== old_dir) {
              errCode = FS.nodePermissions(old_dir, 'w');
              if (errCode) {
                throw new FS.ErrnoError(errCode);
              }
            }
            // remove the node from the lookup hash
            FS.hashRemoveNode(old_node);
            // do the underlying fs rename
            try {
              old_dir.node_ops.rename(old_node, new_dir, new_name);
              // update old node (we do this here to avoid each backend
              // needing to)
              old_node.parent = new_dir;
            } catch (e) {
              throw e;
            } finally {
              // add the node back to the hash (in case node_ops.rename
              // changed its name)
              FS.hashAddNode(old_node);
            }
          },
          rmdir(path) {
            var lookup = FS.lookupPath(path, {parent : true});
            var parent = lookup.node;
            var name = PATH.basename(path);
            var node = FS.lookupNode(parent, name);
            var errCode = FS.mayDelete(parent, name, true);
            if (errCode) {
              throw new FS.ErrnoError(errCode);
            }
            if (!parent.node_ops.rmdir) {
              throw new FS.ErrnoError(63);
            }
            if (FS.isMountpoint(node)) {
              throw new FS.ErrnoError(10);
            }
            parent.node_ops.rmdir(parent, name);
            FS.destroyNode(node);
          },
          readdir(path) {
            var lookup = FS.lookupPath(path, {follow : true});
            var node = lookup.node;
            var readdir = FS.checkOpExists(node.node_ops.readdir, 54);
            return readdir(node);
          },
          unlink(path) {
            var lookup = FS.lookupPath(path, {parent : true});
            var parent = lookup.node;
            if (!parent) {
              throw new FS.ErrnoError(44);
            }
            var name = PATH.basename(path);
            var node = FS.lookupNode(parent, name);
            var errCode = FS.mayDelete(parent, name, false);
            if (errCode) {
              // According to POSIX, we should map EISDIR to EPERM, but
              // we instead do what Linux does (and we must, as we use
              // the musl linux libc).
              throw new FS.ErrnoError(errCode);
            }
            if (!parent.node_ops.unlink) {
              throw new FS.ErrnoError(63);
            }
            if (FS.isMountpoint(node)) {
              throw new FS.ErrnoError(10);
            }
            parent.node_ops.unlink(parent, name);
            FS.destroyNode(node);
          },
          readlink(path) {
            var lookup = FS.lookupPath(path);
            var link = lookup.node;
            if (!link) {
              throw new FS.ErrnoError(44);
            }
            if (!link.node_ops.readlink) {
              throw new FS.ErrnoError(28);
            }
            return link.node_ops.readlink(link);
          },
          stat(path, dontFollow) {
            var lookup = FS.lookupPath(path, {follow : !dontFollow});
            var node = lookup.node;
            var getattr = FS.checkOpExists(node.node_ops.getattr, 63);
            return getattr(node);
          },
          fstat(fd) {
            var stream = FS.getStreamChecked(fd);
            var node = stream.node;
            var getattr = stream.stream_ops.getattr;
            var arg = getattr ? stream : node;
            getattr ??= node.node_ops.getattr;
            FS.checkOpExists(getattr, 63)
            return getattr(arg);
          },
          lstat(path) { return FS.stat(path, true); },
          doChmod(stream, node, mode, dontFollow) {
            FS.doSetAttr(stream, node, {
              mode : (mode & 4095) | (node.mode & ~4095),
              ctime : Date.now(),
              dontFollow
            });
          },
          chmod(path, mode, dontFollow) {
            var node;
            if (typeof path == 'string') {
              var lookup = FS.lookupPath(path, {follow : !dontFollow});
              node = lookup.node;
            } else {
              node = path;
            }
            FS.doChmod(null, node, mode, dontFollow);
          },
          lchmod(path, mode) { FS.chmod(path, mode, true); },
          fchmod(fd, mode) {
            var stream = FS.getStreamChecked(fd);
            FS.doChmod(stream, stream.node, mode, false);
          },
          doChown(stream, node, dontFollow) {
            FS.doSetAttr(stream, node, {
              timestamp : Date.now(),
              dontFollow
              // we ignore the uid / gid for now
            });
          },
          chown(path, uid, gid, dontFollow) {
            var node;
            if (typeof path == 'string') {
              var lookup = FS.lookupPath(path, {follow : !dontFollow});
              node = lookup.node;
            } else {
              node = path;
            }
            FS.doChown(null, node, dontFollow);
          },
          lchown(path, uid, gid) { FS.chown(path, uid, gid, true); },
          fchown(fd, uid, gid) {
            var stream = FS.getStreamChecked(fd);
            FS.doChown(stream, stream.node, false);
          },
          doTruncate(stream, node, len) {
            if (FS.isDir(node.mode)) {
              throw new FS.ErrnoError(31);
            }
            if (!FS.isFile(node.mode)) {
              throw new FS.ErrnoError(28);
            }
            var errCode = FS.nodePermissions(node, 'w');
            if (errCode) {
              throw new FS.ErrnoError(errCode);
            }
            FS.doSetAttr(stream, node, {size : len, timestamp : Date.now()});
          },
          truncate(path, len) {
            if (len < 0) {
              throw new FS.ErrnoError(28);
            }
            var node;
            if (typeof path == 'string') {
              var lookup = FS.lookupPath(path, {follow : true});
              node = lookup.node;
            } else {
              node = path;
            }
            FS.doTruncate(null, node, len);
          },
          ftruncate(fd, len) {
            var stream = FS.getStreamChecked(fd);
            if (len < 0 || (stream.flags & 2097155) === 0) {
              throw new FS.ErrnoError(28);
            }
            FS.doTruncate(stream, stream.node, len);
          },
          utime(path, atime, mtime) {
            var lookup = FS.lookupPath(path, {follow : true});
            var node = lookup.node;
            var setattr = FS.checkOpExists(node.node_ops.setattr, 63);
            setattr(node, {atime : atime, mtime : mtime});
          },
          open(path, flags, mode = 0o666) {
            if (path === "") {
              throw new FS.ErrnoError(44);
            }
            flags =
                typeof flags == 'string' ? FS_modeStringToFlags(flags) : flags;
            if ((flags & 64)) {
              mode = (mode & 4095) | 32768;
            } else {
              mode = 0;
            }
            var node;
            var isDirPath;
            if (typeof path == 'object') {
              node = path;
            } else {
              isDirPath = path.endsWith("/");
              // noent_okay makes it so that if the final component of the path
              // doesn't exist, lookupPath returns `node: undefined`. `path`
              // will be updated to point to the target of all symlinks.
              var lookup = FS.lookupPath(
                  path, {follow : !(flags & 131072), noent_okay : true});
              node = lookup.node;
              path = lookup.path;
            }
            // perhaps we need to create the node
            var created = false;
            if ((flags & 64)) {
              if (node) {
                // if O_CREAT and O_EXCL are set, error out if the node already
                // exists
                if ((flags & 128)) {
                  throw new FS.ErrnoError(20);
                }
              } else if (isDirPath) {
                throw new FS.ErrnoError(31);
              } else {
                // node doesn't exist, try to create it
                // Ignore the permission bits here to ensure we can `open` this
                // new file below. We use chmod below the apply the permissions
                // once the file is open.
                node = FS.mknod(path, mode | 0o777, 0);
                created = true;
              }
            }
            if (!node) {
              throw new FS.ErrnoError(44);
            }
            // can't truncate a device
            if (FS.isChrdev(node.mode)) {
              flags &= ~512;
            }
            // if asked only for a directory, then this must be one
            if ((flags & 65536) && !FS.isDir(node.mode)) {
              throw new FS.ErrnoError(54);
            }
            // check permissions, if this is not a file we just created now (it
            // is ok to create and write to a file with read-only permissions;
            // it is read-only for later use)
            if (!created) {
              var errCode = FS.mayOpen(node, flags);
              if (errCode) {
                throw new FS.ErrnoError(errCode);
              }
            }
            // do truncation if necessary
            if ((flags & 512) && !created) {
              FS.truncate(node, 0);
            }
            // we've already handled these, don't pass down to the underlying
            // vfs
            flags &= ~(128 | 512 | 131072);

            // register the stream with the filesystem
            var stream = FS.createStream({
              node,
              path : FS.getPath(node), // we want the absolute path to the node
              flags,
              seekable : true,
              position : 0,
              stream_ops : node.stream_ops,
              // used by the file family libc calls (fopen, fwrite, ferror,
              // etc.)
              ungotten : [],
              error : false
            });
            // call the new stream's open function
            if (stream.stream_ops.open) {
              stream.stream_ops.open(stream);
            }
            if (created) {
              FS.chmod(node, mode & 0o777);
            }
            if (Module['logReadFiles'] && !(flags & 1)) {
              if (!(path in FS.readFiles)) {
                FS.readFiles[path] = 1;
              }
            }
            return stream;
          },
          close(stream) {
            if (FS.isClosed(stream)) {
              throw new FS.ErrnoError(8);
            }
            if (stream.getdents)
              stream.getdents = null; // free readdir state
            try {
              if (stream.stream_ops.close) {
                stream.stream_ops.close(stream);
              }
            } catch (e) {
              throw e;
            } finally {
              FS.closeStream(stream.fd);
            }
            stream.fd = null;
          },
          isClosed(stream) { return stream.fd === null; },
          llseek(stream, offset, whence) {
            if (FS.isClosed(stream)) {
              throw new FS.ErrnoError(8);
            }
            if (!stream.seekable || !stream.stream_ops.llseek) {
              throw new FS.ErrnoError(70);
            }
            if (whence != 0 && whence != 1 && whence != 2) {
              throw new FS.ErrnoError(28);
            }
            stream.position = stream.stream_ops.llseek(stream, offset, whence);
            stream.ungotten = [];
            return stream.position;
          },
          read(stream, buffer, offset, length, position) {
            assert(offset >= 0);
            if (length < 0 || position < 0) {
              throw new FS.ErrnoError(28);
            }
            if (FS.isClosed(stream)) {
              throw new FS.ErrnoError(8);
            }
            if ((stream.flags & 2097155) === 1) {
              throw new FS.ErrnoError(8);
            }
            if (FS.isDir(stream.node.mode)) {
              throw new FS.ErrnoError(31);
            }
            if (!stream.stream_ops.read) {
              throw new FS.ErrnoError(28);
            }
            var seeking = typeof position != 'undefined';
            if (!seeking) {
              position = stream.position;
            } else if (!stream.seekable) {
              throw new FS.ErrnoError(70);
            }
            var bytesRead = stream.stream_ops.read(stream, buffer, offset,
                                                   length, position);
            if (!seeking)
              stream.position += bytesRead;
            return bytesRead;
          },
          write(stream, buffer, offset, length, position, canOwn) {
            assert(offset >= 0);
            if (length < 0 || position < 0) {
              throw new FS.ErrnoError(28);
            }
            if (FS.isClosed(stream)) {
              throw new FS.ErrnoError(8);
            }
            if ((stream.flags & 2097155) === 0) {
              throw new FS.ErrnoError(8);
            }
            if (FS.isDir(stream.node.mode)) {
              throw new FS.ErrnoError(31);
            }
            if (!stream.stream_ops.write) {
              throw new FS.ErrnoError(28);
            }
            if (stream.seekable && stream.flags & 1024) {
              // seek to the end before writing in append mode
              FS.llseek(stream, 0, 2);
            }
            var seeking = typeof position != 'undefined';
            if (!seeking) {
              position = stream.position;
            } else if (!stream.seekable) {
              throw new FS.ErrnoError(70);
            }
            var bytesWritten = stream.stream_ops.write(
                stream, buffer, offset, length, position, canOwn);
            if (!seeking)
              stream.position += bytesWritten;
            return bytesWritten;
          },
          mmap(stream, length, position, prot, flags) {
            // User requests writing to file (prot & PROT_WRITE != 0).
            // Checking if we have permissions to write to the file unless
            // MAP_PRIVATE flag is set. According to POSIX spec it is possible
            // to write to file opened in read-only mode with MAP_PRIVATE flag,
            // as all modifications will be visible only in the memory of
            // the current process.
            if ((prot & 2) !== 0 && (flags & 2) === 0 &&
                (stream.flags & 2097155) !== 2) {
              throw new FS.ErrnoError(2);
            }
            if ((stream.flags & 2097155) === 1) {
              throw new FS.ErrnoError(2);
            }
            if (!stream.stream_ops.mmap) {
              throw new FS.ErrnoError(43);
            }
            if (!length) {
              throw new FS.ErrnoError(28);
            }
            return stream.stream_ops.mmap(stream, length, position, prot,
                                          flags);
          },
          msync(stream, buffer, offset, length, mmapFlags) {
            assert(offset >= 0);
            if (!stream.stream_ops.msync) {
              return 0;
            }
            return stream.stream_ops.msync(stream, buffer, offset, length,
                                           mmapFlags);
          },
          ioctl(stream, cmd, arg) {
            if (!stream.stream_ops.ioctl) {
              throw new FS.ErrnoError(59);
            }
            return stream.stream_ops.ioctl(stream, cmd, arg);
          },
          readFile(path, opts = {}) {
            opts.flags = opts.flags || 0;
            opts.encoding = opts.encoding || 'binary';
            if (opts.encoding !== 'utf8' && opts.encoding !== 'binary') {
              throw new Error(`Invalid encoding type "${opts.encoding}"`);
            }
            var stream = FS.open(path, opts.flags);
            var stat = FS.stat(path);
            var length = stat.size;
            var buf = new Uint8Array(length);
            FS.read(stream, buf, 0, length, 0);
            if (opts.encoding === 'utf8') {
              buf = UTF8ArrayToString(buf);
            }
            FS.close(stream);
            return buf;
          },
          writeFile(path, data, opts = {}) {
            opts.flags = opts.flags || 577;
            var stream = FS.open(path, opts.flags, opts.mode);
            if (typeof data == 'string') {
              data = new Uint8Array(intArrayFromString(data, true));
            }
            if (ArrayBuffer.isView(data)) {
              FS.write(stream, data, 0, data.byteLength, undefined,
                       opts.canOwn);
            } else {
              throw new Error('Unsupported data type');
            }
            FS.close(stream);
          },
          cwd : () => FS.currentPath,
          chdir(path) {
            var lookup = FS.lookupPath(path, {follow : true});
            if (lookup.node === null) {
              throw new FS.ErrnoError(44);
            }
            if (!FS.isDir(lookup.node.mode)) {
              throw new FS.ErrnoError(54);
            }
            var errCode = FS.nodePermissions(lookup.node, 'x');
            if (errCode) {
              throw new FS.ErrnoError(errCode);
            }
            FS.currentPath = lookup.path;
          },
          createDefaultDirectories() {
            FS.mkdir('/tmp');
            FS.mkdir('/home');
            FS.mkdir('/home/web_user');
          },
          createDefaultDevices() {
            // create /dev
            FS.mkdir('/dev');
            // setup /dev/null
            FS.registerDevice(FS.makedev(1, 3), {
              read : () => 0,
              write : (stream, buffer, offset, length, pos) => length,
              llseek : () => 0,
            });
            FS.mkdev('/dev/null', FS.makedev(1, 3));
            // setup /dev/tty and /dev/tty1
            // stderr needs to print output using err() rather than out()
            // so we register a second tty just for it.
            TTY.register(FS.makedev(5, 0), TTY.default_tty_ops);
            TTY.register(FS.makedev(6, 0), TTY.default_tty1_ops);
            FS.mkdev('/dev/tty', FS.makedev(5, 0));
            FS.mkdev('/dev/tty1', FS.makedev(6, 0));
            // setup /dev/[u]random
            // use a buffer to avoid overhead of individual crypto calls per
            // byte
            var randomBuffer = new Uint8Array(1024), randomLeft = 0;
            var randomByte = () => {
              if (randomLeft === 0) {
                randomFill(randomBuffer);
                randomLeft = randomBuffer.byteLength;
              }
              return randomBuffer[--randomLeft];
            };
            FS.createDevice('/dev', 'random', randomByte);
            FS.createDevice('/dev', 'urandom', randomByte);
            // we're not going to emulate the actual shm device,
            // just create the tmp dirs that reside in it commonly
            FS.mkdir('/dev/shm');
            FS.mkdir('/dev/shm/tmp');
          },
          createSpecialDirectories() {
            // create /proc/self/fd which allows /proc/self/fd/6 => readlink
            // gives the name of the stream for fd 6 (see test_unistd_ttyname)
            FS.mkdir('/proc');
            var proc_self = FS.mkdir('/proc/self');
            FS.mkdir('/proc/self/fd');
            FS.mount({
              mount() {
                var node = FS.createNode(proc_self, 'fd', 16895, 73);
                node.stream_ops = {
                  llseek : MEMFS.stream_ops.llseek,
                };
                node.node_ops = {
                  lookup(parent, name) {
                    var fd = +name;
                    var stream = FS.getStreamChecked(fd);
                    var ret = {
                      parent : null,
                      mount : {mountpoint : 'fake'},
                      node_ops : {readlink : () => stream.path},
                      id : fd + 1,
                    };
                    ret.parent = ret; // make it look like a simple root node
                    return ret;
                  },
                  readdir() {
                    return Array.from(FS.streams.entries())
                        .filter(([ k, v ]) => v)
                        .map(([ k, v ]) => k.toString());
                  }
                };
                return node;
              }
            },
                     {}, '/proc/self/fd');
          },
          createStandardStreams(input, output, error) {
            // TODO deprecate the old functionality of a single
            // input / output callback and that utilizes FS.createDevice
            // and instead require a unique set of stream ops

            // by default, we symlink the standard streams to the
            // default tty devices. however, if the standard streams
            // have been overwritten we create a unique device for
            // them instead.
            if (input) {
              FS.createDevice('/dev', 'stdin', input);
            } else {
              FS.symlink('/dev/tty', '/dev/stdin');
            }
            if (output) {
              FS.createDevice('/dev', 'stdout', null, output);
            } else {
              FS.symlink('/dev/tty', '/dev/stdout');
            }
            if (error) {
              FS.createDevice('/dev', 'stderr', null, error);
            } else {
              FS.symlink('/dev/tty1', '/dev/stderr');
            }

            // open default streams for the stdin, stdout and stderr devices
            var stdin = FS.open('/dev/stdin', 0);
            var stdout = FS.open('/dev/stdout', 1);
            var stderr = FS.open('/dev/stderr', 1);
            assert(stdin.fd === 0, `invalid handle for stdin (${stdin.fd})`);
            assert(stdout.fd === 1, `invalid handle for stdout (${stdout.fd})`);
            assert(stderr.fd === 2, `invalid handle for stderr (${stderr.fd})`);
          },
          staticInit() {
            FS.nameTable = new Array(4096);

            FS.mount(MEMFS, {}, '/');

            FS.createDefaultDirectories();
            FS.createDefaultDevices();
            FS.createSpecialDirectories();

            FS.filesystems = {
              'MEMFS' : MEMFS,
              'WORKERFS' : WORKERFS,
              'PROXYFS' : PROXYFS,
            };
          },
          init(input, output, error) {
            assert(
                !FS.initialized,
                'FS.init was previously called. If you want to initialize later with custom parameters, remove any earlier calls (note that one is automatically added to the generated code)');
            FS.initialized = true;

            // Allow Module.stdin etc. to provide defaults, if none explicitly
            // passed to us here
            input ??= Module['stdin'];
            output ??= Module['stdout'];
            error ??= Module['stderr'];

            FS.createStandardStreams(input, output, error);
          },
          quit() {
            FS.initialized = false;
            // force-flush all streams, so we get musl std streams printed out
            _fflush(0);
            // close all of our streams
            for (var stream of FS.streams) {
              if (stream) {
                FS.close(stream);
              }
            }
          },
          findObject(path, dontResolveLastLink) {
            var ret = FS.analyzePath(path, dontResolveLastLink);
            if (!ret.exists) {
              return null;
            }
            return ret.object;
          },
          analyzePath(path, dontResolveLastLink) {
            // operate from within the context of the symlink's target
            try {
              var lookup = FS.lookupPath(path, {follow : !dontResolveLastLink});
              path = lookup.path;
            } catch (e) {
            }
            var ret = {
              isRoot : false,
              exists : false,
              error : 0,
              name : null,
              path : null,
              object : null,
              parentExists : false,
              parentPath : null,
              parentObject : null
            };
            try {
              var lookup = FS.lookupPath(path, {parent : true});
              ret.parentExists = true;
              ret.parentPath = lookup.path;
              ret.parentObject = lookup.node;
              ret.name = PATH.basename(path);
              lookup = FS.lookupPath(path, {follow : !dontResolveLastLink});
              ret.exists = true;
              ret.path = lookup.path;
              ret.object = lookup.node;
              ret.name = lookup.node.name;
              ret.isRoot = lookup.path === '/';
            } catch (e) {
              ret.error = e.errno;
            };
            return ret;
          },
          createPath(parent, path, canRead, canWrite) {
            parent = typeof parent == 'string' ? parent : FS.getPath(parent);
            var parts = path.split('/').reverse();
            while (parts.length) {
              var part = parts.pop();
              if (!part)
                continue;
              var current = PATH.join2(parent, part);
              try {
                FS.mkdir(current);
              } catch (e) {
                if (e.errno != 20)
                  throw e;
              }
              parent = current;
            }
            return current;
          },
          createFile(parent, name, properties, canRead, canWrite) {
            var path = PATH.join2(
                typeof parent == 'string' ? parent : FS.getPath(parent), name);
            var mode = FS_getMode(canRead, canWrite);
            return FS.create(path, mode);
          },
          createDataFile(parent, name, data, canRead, canWrite, canOwn) {
            var path = name;
            if (parent) {
              parent = typeof parent == 'string' ? parent : FS.getPath(parent);
              path = name ? PATH.join2(parent, name) : parent;
            }
            var mode = FS_getMode(canRead, canWrite);
            var node = FS.create(path, mode);
            if (data) {
              if (typeof data == 'string') {
                var arr = new Array(data.length);
                for (var i = 0, len = data.length; i < len; ++i)
                  arr[i] = data.charCodeAt(i);
                data = arr;
              }
              // make sure we can write to the file
              FS.chmod(node, mode | 146);
              var stream = FS.open(node, 577);
              FS.write(stream, data, 0, data.length, 0, canOwn);
              FS.close(stream);
              FS.chmod(node, mode);
            }
          },
          createDevice(parent, name, input, output) {
            var path = PATH.join2(
                typeof parent == 'string' ? parent : FS.getPath(parent), name);
            var mode = FS_getMode(!!input, !!output);
            FS.createDevice.major ??= 64;
            var dev = FS.makedev(FS.createDevice.major++, 0);
            // Create a fake device that a set of stream ops to emulate
            // the old behavior.
            FS.registerDevice(dev, {
              open(stream) { stream.seekable = false; },
              close(stream) {
                // flush any pending line data
                if (output?.buffer?.length) {
                  output(10);
                }
              },
              read(stream, buffer, offset, length, pos /* ignored */) {
                var bytesRead = 0;
                for (var i = 0; i < length; i++) {
                  var result;
                  try {
                    result = input();
                  } catch (e) {
                    throw new FS.ErrnoError(29);
                  }
                  if (result === undefined && bytesRead === 0) {
                    throw new FS.ErrnoError(6);
                  }
                  if (result === null || result === undefined)
                    break;
                  bytesRead++;
                  buffer[offset + i] = result;
                }
                if (bytesRead) {
                  stream.node.atime = Date.now();
                }
                return bytesRead;
              },
              write(stream, buffer, offset, length, pos) {
                for (var i = 0; i < length; i++) {
                  try {
                    output(buffer[offset + i]);
                  } catch (e) {
                    throw new FS.ErrnoError(29);
                  }
                }
                if (length) {
                  stream.node.mtime = stream.node.ctime = Date.now();
                }
                return i;
              }
            });
            return FS.mkdev(path, mode, dev);
          },
          forceLoadFile(obj) {
            if (obj.isDevice || obj.isFolder || obj.link || obj.contents)
              return true;
            if (typeof XMLHttpRequest != 'undefined') {
              throw new Error(
                  "Lazy loading should have been performed (contents set) in createLazyFile, but it was not. Lazy loading only works in web workers. Use --embed-file or --preload-file in emcc on the main thread.");
            } else { // Command-line.
              try {
                obj.contents = readBinary(obj.url);
                obj.usedBytes = obj.contents.length;
              } catch (e) {
                throw new FS.ErrnoError(29);
              }
            }
          },
          createLazyFile(parent, name, url, canRead, canWrite) {
            // Lazy chunked Uint8Array (implements get and length from
            // Uint8Array). Actual getting is abstracted away for eventual
            // reuse.
            class LazyUint8Array {
              lengthKnown = false;
              chunks = []; // Loaded chunks. Index is the chunk number
              get(idx) {
                if (idx > this.length - 1 || idx < 0) {
                  return undefined;
                }
                var chunkOffset = idx % this.chunkSize;
                var chunkNum = (idx / this.chunkSize) | 0;
                return this.getter(chunkNum)[chunkOffset];
              }
              setDataGetter(getter) { this.getter = getter; }
              cacheLength() {
                // Find length
                var xhr = new XMLHttpRequest();
                xhr.open('HEAD', url, false);
                xhr.send(null);
                if (!(xhr.status >= 200 && xhr.status < 300 ||
                      xhr.status === 304))
                  throw new Error("Couldn't load " + url +
                                  ". Status: " + xhr.status);
                var datalength =
                    Number(xhr.getResponseHeader("Content-length"));
                var header;
                var hasByteServing =
                    (header = xhr.getResponseHeader("Accept-Ranges")) &&
                    header === "bytes";
                var usesGzip =
                    (header = xhr.getResponseHeader("Content-Encoding")) &&
                    header === "gzip";

                var chunkSize = 1024 * 1024; // Chunk size in bytes

                if (!hasByteServing)
                  chunkSize = datalength;

                // Function to get a range from the remote URL.
                var doXHR = (from, to) => {
                  if (from > to)
                    throw new Error("invalid range (" + from + ", " + to +
                                    ") or no bytes requested!");
                  if (to > datalength - 1)
                    throw new Error("only " + datalength +
                                    " bytes available! programmer error!");

                  // TODO: Use mozResponseArrayBuffer, responseStream, etc. if
                  // available.
                  var xhr = new XMLHttpRequest();
                  xhr.open('GET', url, false);
                  if (datalength !== chunkSize)
                    xhr.setRequestHeader("Range", "bytes=" + from + "-" + to);

                  // Some hints to the browser that we want binary data.
                  xhr.responseType = 'arraybuffer';
                  if (xhr.overrideMimeType) {
                    xhr.overrideMimeType('text/plain; charset=x-user-defined');
                  }

                  xhr.send(null);
                  if (!(xhr.status >= 200 && xhr.status < 300 ||
                        xhr.status === 304))
                    throw new Error("Couldn't load " + url +
                                    ". Status: " + xhr.status);
                  if (xhr.response !== undefined) {
                    return new Uint8Array(
                        /** @type{Array<number>} */ (xhr.response || []));
                  }
                  return intArrayFromString(xhr.responseText || '', true);
                };
                var lazyArray = this;
                lazyArray.setDataGetter((chunkNum) => {
                  var start = chunkNum * chunkSize;
                  var end =
                      (chunkNum + 1) * chunkSize - 1; // including this byte
                  end = Math.min(end, datalength -
                                          1); // if datalength-1 is selected,
                                              // this is the last block
                  if (typeof lazyArray.chunks[chunkNum] == 'undefined') {
                    lazyArray.chunks[chunkNum] = doXHR(start, end);
                  }
                  if (typeof lazyArray.chunks[chunkNum] == 'undefined')
                    throw new Error('doXHR failed!');
                  return lazyArray.chunks[chunkNum];
                });

                if (usesGzip || !datalength) {
                  // if the server uses gzip or doesn't supply the length, we
                  // have to download the whole file to get the (uncompressed)
                  // length
                  chunkSize = datalength = 1; // this will force getter(0)/doXHR
                                              // do download the whole file
                  datalength = this.getter(0).length;
                  chunkSize = datalength;
                  out("LazyFiles on gzip forces download of the whole file when length is accessed");
                }

                this._length = datalength;
                this._chunkSize = chunkSize;
                this.lengthKnown = true;
              }
              get length() {
                if (!this.lengthKnown) {
                  this.cacheLength();
                }
                return this._length;
              }
              get chunkSize() {
                if (!this.lengthKnown) {
                  this.cacheLength();
                }
                return this._chunkSize;
              }
            }

            if (typeof XMLHttpRequest != 'undefined') {
              if (!ENVIRONMENT_IS_WORKER)
                throw 'Cannot do synchronous binary XHRs outside webworkers in modern browsers. Use --embed-file or --preload-file in emcc';
              var lazyArray = new LazyUint8Array();
              var properties = {isDevice : false, contents : lazyArray};
            } else {
              var properties = {isDevice : false, url : url};
            }

            var node =
                FS.createFile(parent, name, properties, canRead, canWrite);
            // This is a total hack, but I want to get this lazy file code out
            // of the core of MEMFS. If we want to keep this lazy file concept I
            // feel it should be its own thin LAZYFS proxying calls to MEMFS.
            if (properties.contents) {
              node.contents = properties.contents;
            } else if (properties.url) {
              node.contents = null;
              node.url = properties.url;
            }
            // Add a function that defers querying the file size until it is
            // asked the first time.
            Object.defineProperties(node, {
              usedBytes : {get : function() { return this.contents.length; }}
            });
            // override each stream op with one that tries to force load the
            // lazy file first
            var stream_ops = {};
            var keys = Object.keys(node.stream_ops);
            keys.forEach((key) => {
              var fn = node.stream_ops[key];
              stream_ops[key] = (...args) => {
                FS.forceLoadFile(node);
                return fn(...args);
              };
            });
            function writeChunks(stream, buffer, offset, length, position) {
              var contents = stream.node.contents;
              if (position >= contents.length)
                return 0;
              var size = Math.min(contents.length - position, length);
              assert(size >= 0);
              if (contents.slice) { // normal array
                for (var i = 0; i < size; i++) {
                  buffer[offset + i] = contents[position + i];
                }
              } else {
                for (var i = 0; i < size;
                     i++) { // LazyUint8Array from sync binary XHR
                  buffer[offset + i] = contents.get(position + i);
                }
              }
              return size;
            }
            // use a custom read function
            stream_ops.read = (stream, buffer, offset, length, position) => {
              FS.forceLoadFile(node);
              return writeChunks(stream, buffer, offset, length, position)
            };
            // use a custom mmap function
            stream_ops.mmap = (stream, length, position, prot, flags) => {
              FS.forceLoadFile(node);
              var ptr = mmapAlloc(length);
              if (!ptr) {
                throw new FS.ErrnoError(48);
              }
              writeChunks(stream, HEAP8, ptr, length, position);
              return {ptr, allocated : true};
            };
            node.stream_ops = stream_ops;
            return node;
          },
          absolutePath() {
            abort(
                'FS.absolutePath has been removed; use PATH_FS.resolve instead');
          },
          createFolder() {
            abort('FS.createFolder has been removed; use FS.mkdir instead');
          },
          createLink() {
            abort('FS.createLink has been removed; use FS.symlink instead');
          },
          joinPath() {
            abort('FS.joinPath has been removed; use PATH.join instead');
          },
          mmapAlloc() {
            abort(
                'FS.mmapAlloc has been replaced by the top level function mmapAlloc');
          },
          standardizePath() {
            abort(
                'FS.standardizePath has been removed; use PATH.normalize instead');
          },
    };

    var SYSCALLS = {
      DEFAULT_POLLMASK : 5,
      calculateAt(dirfd, path, allowEmpty) {
        if (PATH.isAbs(path)) {
          return path;
        }
        // relative path
        var dir;
        if (dirfd === -100) {
          dir = FS.cwd();
        } else {
          var dirstream = SYSCALLS.getStreamFromFD(dirfd);
          dir = dirstream.path;
        }
        if (path.length == 0) {
          if (!allowEmpty) {
            throw new FS.ErrnoError(44);
            ;
          }
          return dir;
        }
        return dir + '/' + path;
      },
      writeStat(buf, stat) {
        HEAPU32[((buf) >> 2)] = stat.dev;
        HEAPU32[(((buf) + (4)) >> 2)] = stat.mode;
        HEAPU32[(((buf) + (8)) >> 2)] = stat.nlink;
        HEAPU32[(((buf) + (12)) >> 2)] = stat.uid;
        HEAPU32[(((buf) + (16)) >> 2)] = stat.gid;
        HEAPU32[(((buf) + (20)) >> 2)] = stat.rdev;
        HEAP64[(((buf) + (24)) >> 3)] = BigInt(stat.size);
        HEAP32[(((buf) + (32)) >> 2)] = 4096;
        HEAP32[(((buf) + (36)) >> 2)] = stat.blocks;
        var atime = stat.atime.getTime();
        var mtime = stat.mtime.getTime();
        var ctime = stat.ctime.getTime();
        HEAP64[(((buf) + (40)) >> 3)] = BigInt(Math.floor(atime / 1000));
        HEAPU32[(((buf) + (48)) >> 2)] = (atime % 1000) * 1000 * 1000;
        HEAP64[(((buf) + (56)) >> 3)] = BigInt(Math.floor(mtime / 1000));
        HEAPU32[(((buf) + (64)) >> 2)] = (mtime % 1000) * 1000 * 1000;
        HEAP64[(((buf) + (72)) >> 3)] = BigInt(Math.floor(ctime / 1000));
        HEAPU32[(((buf) + (80)) >> 2)] = (ctime % 1000) * 1000 * 1000;
        HEAP64[(((buf) + (88)) >> 3)] = BigInt(stat.ino);
        return 0;
      },
      writeStatFs(buf, stats) {
        HEAPU32[(((buf) + (4)) >> 2)] = stats.bsize;
        HEAPU32[(((buf) + (60)) >> 2)] = stats.bsize;
        HEAP64[(((buf) + (8)) >> 3)] = BigInt(stats.blocks);
        HEAP64[(((buf) + (16)) >> 3)] = BigInt(stats.bfree);
        HEAP64[(((buf) + (24)) >> 3)] = BigInt(stats.bavail);
        HEAP64[(((buf) + (32)) >> 3)] = BigInt(stats.files);
        HEAP64[(((buf) + (40)) >> 3)] = BigInt(stats.ffree);
        HEAPU32[(((buf) + (48)) >> 2)] = stats.fsid;
        HEAPU32[(((buf) + (64)) >> 2)] = stats.flags; // ST_NOSUID
        HEAPU32[(((buf) + (56)) >> 2)] = stats.namelen;
      },
      doMsync(addr, stream, len, flags, offset) {
        if (!FS.isFile(stream.node.mode)) {
          throw new FS.ErrnoError(43);
        }
        if (flags & 2) {
          // MAP_PRIVATE calls need not to be synced back to underlying fs
          return 0;
        }
        var buffer = HEAPU8.slice(addr, addr + len);
        FS.msync(stream, buffer, offset, len, flags);
      },
      getStreamFromFD(fd) {
        var stream = FS.getStreamChecked(fd);
        return stream;
      },
      varargs : undefined,
      getStr(ptr) {
        var ret = UTF8ToString(ptr);
        return ret;
      },
    };
    function ___syscall_fcntl64(fd, cmd, varargs) {
      SYSCALLS.varargs = varargs;
      try {

        var stream = SYSCALLS.getStreamFromFD(fd);
        switch (cmd) {
        case 0: {
          var arg = syscallGetVarargI();
          if (arg < 0) {
            return -28;
          }
          while (FS.streams[arg]) {
            arg++;
          }
          var newStream;
          newStream = FS.dupStream(stream, arg);
          return newStream.fd;
        }
        case 1:
        case 2:
          return 0; // FD_CLOEXEC makes no sense for a single process.
        case 3:
          return stream.flags;
        case 4: {
          var arg = syscallGetVarargI();
          stream.flags |= arg;
          return 0;
        }
        case 12: {
          var arg = syscallGetVarargP();
          var offset = 0;
          // We're always unlocked.
          HEAP16[(((arg) + (offset)) >> 1)] = 2;
          return 0;
        }
        case 13:
        case 14:
          // Pretend that the locking is successful. These are process-level
          // locks, and Emscripten programs are a single process. If we
          // supported linking a filesystem between programs, we'd need to do
          // more here. See
          // https://github.com/emscripten-core/emscripten/issues/23697
          return 0;
        }
        return -28;
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return -e.errno;
      }
    }

    function ___syscall_fstat64(fd, buf) {
      try {

        return SYSCALLS.writeStat(buf, FS.fstat(fd));
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return -e.errno;
      }
    }

    var stringToUTF8 = (str, outPtr, maxBytesToWrite) => {
      assert(
          typeof maxBytesToWrite == 'number',
          'stringToUTF8(str, outPtr, maxBytesToWrite) is missing the third parameter that specifies the length of the output buffer!');
      return stringToUTF8Array(str, HEAPU8, outPtr, maxBytesToWrite);
    };
    function ___syscall_getcwd(buf, size) {
      try {

        if (size === 0)
          return -28;
        var cwd = FS.cwd();
        var cwdLengthInBytes = lengthBytesUTF8(cwd) + 1;
        if (size < cwdLengthInBytes)
          return -68;
        stringToUTF8(cwd, buf, size);
        return cwdLengthInBytes;
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return -e.errno;
      }
    }

    function ___syscall_ioctl(fd, op, varargs) {
      SYSCALLS.varargs = varargs;
      try {

        var stream = SYSCALLS.getStreamFromFD(fd);
        switch (op) {
        case 21509: {
          if (!stream.tty)
            return -59;
          return 0;
        }
        case 21505: {
          if (!stream.tty)
            return -59;
          if (stream.tty.ops.ioctl_tcgets) {
            var termios = stream.tty.ops.ioctl_tcgets(stream);
            var argp = syscallGetVarargP();
            HEAP32[((argp) >> 2)] = termios.c_iflag || 0;
            HEAP32[(((argp) + (4)) >> 2)] = termios.c_oflag || 0;
            HEAP32[(((argp) + (8)) >> 2)] = termios.c_cflag || 0;
            HEAP32[(((argp) + (12)) >> 2)] = termios.c_lflag || 0;
            for (var i = 0; i < 32; i++) {
              HEAP8[(argp + i) + (17)] = termios.c_cc[i] || 0;
            }
            return 0;
          }
          return 0;
        }
        case 21510:
        case 21511:
        case 21512: {
          if (!stream.tty)
            return -59;
          return 0; // no-op, not actually adjusting terminal settings
        }
        case 21506:
        case 21507:
        case 21508: {
          if (!stream.tty)
            return -59;
          if (stream.tty.ops.ioctl_tcsets) {
            var argp = syscallGetVarargP();
            var c_iflag = HEAP32[((argp) >> 2)];
            var c_oflag = HEAP32[(((argp) + (4)) >> 2)];
            var c_cflag = HEAP32[(((argp) + (8)) >> 2)];
            var c_lflag = HEAP32[(((argp) + (12)) >> 2)];
            var c_cc = [] for (var i = 0; i < 32; i++) {
              c_cc.push(HEAP8[(argp + i) + (17)]);
            }
            return stream.tty.ops.ioctl_tcsets(
                stream.tty, op, {c_iflag, c_oflag, c_cflag, c_lflag, c_cc});
          }
          return 0; // no-op, not actually adjusting terminal settings
        }
        case 21519: {
          if (!stream.tty)
            return -59;
          var argp = syscallGetVarargP();
          HEAP32[((argp) >> 2)] = 0;
          return 0;
        }
        case 21520: {
          if (!stream.tty)
            return -59;
          return -28; // not supported
        }
        case 21537:
        case 21531: {
          var argp = syscallGetVarargP();
          return FS.ioctl(stream, op, argp);
        }
        case 21523: {
          // TODO: in theory we should write to the winsize struct that gets
          // passed in, but for now musl doesn't read anything on it
          if (!stream.tty)
            return -59;
          if (stream.tty.ops.ioctl_tiocgwinsz) {
            var winsize = stream.tty.ops.ioctl_tiocgwinsz(stream.tty);
            var argp = syscallGetVarargP();
            HEAP16[((argp) >> 1)] = winsize[0];
            HEAP16[(((argp) + (2)) >> 1)] = winsize[1];
          }
          return 0;
        }
        case 21524: {
          // TODO: technically, this ioctl call should change the window size.
          // but, since emscripten doesn't have any concept of a terminal window
          // yet, we'll just silently throw it away as we do TIOCGWINSZ
          if (!stream.tty)
            return -59;
          return 0;
        }
        case 21515: {
          if (!stream.tty)
            return -59;
          return 0;
        }
        default:
          return -28; // not supported
        }
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return -e.errno;
      }
    }

    function ___syscall_lstat64(path, buf) {
      try {

        path = SYSCALLS.getStr(path);
        return SYSCALLS.writeStat(buf, FS.lstat(path));
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return -e.errno;
      }
    }

    function ___syscall_newfstatat(dirfd, path, buf, flags) {
      try {

        path = SYSCALLS.getStr(path);
        var nofollow = flags & 256;
        var allowEmpty = flags & 4096;
        flags = flags & (~6400);
        assert(!flags, `unknown flags in __syscall_newfstatat: ${flags}`);
        path = SYSCALLS.calculateAt(dirfd, path, allowEmpty);
        return SYSCALLS.writeStat(buf,
                                  nofollow ? FS.lstat(path) : FS.stat(path));
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return -e.errno;
      }
    }

    function ___syscall_openat(dirfd, path, flags, varargs) {
      SYSCALLS.varargs = varargs;
      try {

        path = SYSCALLS.getStr(path);
        path = SYSCALLS.calculateAt(dirfd, path);
        var mode = varargs ? syscallGetVarargI() : 0;
        return FS.open(path, flags, mode).fd;
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return -e.errno;
      }
    }

    function ___syscall_rmdir(path) {
      try {

        path = SYSCALLS.getStr(path);
        FS.rmdir(path);
        return 0;
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return -e.errno;
      }
    }

    function ___syscall_stat64(path, buf) {
      try {

        path = SYSCALLS.getStr(path);
        return SYSCALLS.writeStat(buf, FS.stat(path));
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return -e.errno;
      }
    }

    function ___syscall_unlinkat(dirfd, path, flags) {
      try {

        path = SYSCALLS.getStr(path);
        path = SYSCALLS.calculateAt(dirfd, path);
        if (!flags) {
          FS.unlink(path);
        } else if (flags === 512) {
          FS.rmdir(path);
        } else {
          return -28;
        }
        return 0;
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return -e.errno;
      }
    }

    var getCppExceptionTag = () => ___cpp_exception;

    var getCppExceptionThrownObjectFromWebAssemblyException = (ex) => {
      // In Wasm EH, the value extracted from WebAssembly.Exception is a pointer
      // to the unwind header. Convert it to the actual thrown value.
      var unwind_header = ex.getArg(getCppExceptionTag(), 0);
      return ___thrown_object_from_unwind_exception(unwind_header);
    };

    var stackSave = () => _emscripten_stack_get_current();

    var stackRestore = (val) => __emscripten_stack_restore(val);

    var stackAlloc = (sz) => __emscripten_stack_alloc(sz);

    var getExceptionMessageCommon = (ptr) => {
      var sp = stackSave();
      var type_addr_addr = stackAlloc(4);
      var message_addr_addr = stackAlloc(4);
      ___get_exception_message(ptr, type_addr_addr, message_addr_addr);
      var type_addr = HEAPU32[((type_addr_addr) >> 2)];
      var message_addr = HEAPU32[((message_addr_addr) >> 2)];
      var type = UTF8ToString(type_addr);
      _free(type_addr);
      var message;
      if (message_addr) {
        message = UTF8ToString(message_addr);
        _free(message_addr);
      }
      stackRestore(sp);
      return [ type, message ];
    };
    var getExceptionMessage = (ex) => {
      var ptr = getCppExceptionThrownObjectFromWebAssemblyException(ex);
      return getExceptionMessageCommon(ptr);
    };
    var ___throw_exception_with_stack_trace = (ex) => {
      var e = new WebAssembly.Exception(getCppExceptionTag(), [ ex ],
                                        {traceStack : true});
      e.message = getExceptionMessage(e);
      throw e;
    };

    var __abort_js = () => abort('native code called abort()');

    var runtimeKeepaliveCounter = 0;
    var __emscripten_runtime_keepalive_clear = () => {
      noExitRuntime = false;
      runtimeKeepaliveCounter = 0;
    };

    var INT53_MAX = 9007199254740992;

    var INT53_MIN = -9007199254740992;
    var bigintToI53Checked = (num) =>
        (num < INT53_MIN || num > INT53_MAX) ? NaN : Number(num);
    function __gmtime_js(time, tmPtr) {
      time = bigintToI53Checked(time);

      var date = new Date(time * 1000);
      HEAP32[((tmPtr) >> 2)] = date.getUTCSeconds();
      HEAP32[(((tmPtr) + (4)) >> 2)] = date.getUTCMinutes();
      HEAP32[(((tmPtr) + (8)) >> 2)] = date.getUTCHours();
      HEAP32[(((tmPtr) + (12)) >> 2)] = date.getUTCDate();
      HEAP32[(((tmPtr) + (16)) >> 2)] = date.getUTCMonth();
      HEAP32[(((tmPtr) + (20)) >> 2)] = date.getUTCFullYear() - 1900;
      HEAP32[(((tmPtr) + (24)) >> 2)] = date.getUTCDay();
      var start = Date.UTC(date.getUTCFullYear(), 0, 1, 0, 0, 0, 0);
      var yday = ((date.getTime() - start) / (1000 * 60 * 60 * 24)) | 0;
      HEAP32[(((tmPtr) + (28)) >> 2)] = yday;
      ;
    }

    var isLeapYear = (year) =>
        year % 4 === 0 && (year % 100 !== 0 || year % 400 === 0);

    var MONTH_DAYS_LEAP_CUMULATIVE =
        [ 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 ];

    var MONTH_DAYS_REGULAR_CUMULATIVE =
        [ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 ];
    var ydayFromDate = (date) => {
      var leap = isLeapYear(date.getFullYear());
      var monthDaysCumulative =
          (leap ? MONTH_DAYS_LEAP_CUMULATIVE : MONTH_DAYS_REGULAR_CUMULATIVE);
      var yday = monthDaysCumulative[date.getMonth()] + date.getDate() -
                 1; // -1 since it's days since Jan 1

      return yday;
    };

    function __localtime_js(time, tmPtr) {
      time = bigintToI53Checked(time);

      var date = new Date(time * 1000);
      HEAP32[((tmPtr) >> 2)] = date.getSeconds();
      HEAP32[(((tmPtr) + (4)) >> 2)] = date.getMinutes();
      HEAP32[(((tmPtr) + (8)) >> 2)] = date.getHours();
      HEAP32[(((tmPtr) + (12)) >> 2)] = date.getDate();
      HEAP32[(((tmPtr) + (16)) >> 2)] = date.getMonth();
      HEAP32[(((tmPtr) + (20)) >> 2)] = date.getFullYear() - 1900;
      HEAP32[(((tmPtr) + (24)) >> 2)] = date.getDay();

      var yday = ydayFromDate(date) | 0;
      HEAP32[(((tmPtr) + (28)) >> 2)] = yday;
      HEAP32[(((tmPtr) + (36)) >> 2)] = -(date.getTimezoneOffset() * 60);

      // Attention: DST is in December in South, and some regions don't have DST
      // at all.
      var start = new Date(date.getFullYear(), 0, 1);
      var summerOffset = new Date(date.getFullYear(), 6, 1).getTimezoneOffset();
      var winterOffset = start.getTimezoneOffset();
      var dst =
          (summerOffset != winterOffset &&
           date.getTimezoneOffset() == Math.min(winterOffset, summerOffset)) |
          0;
      HEAP32[(((tmPtr) + (32)) >> 2)] = dst;
      ;
    }

    var __tzset_js = (timezone, daylight, std_name, dst_name) => {
      // TODO: Use (malleable) environment variables instead of system settings.
      var currentYear = new Date().getFullYear();
      var winter = new Date(currentYear, 0, 1);
      var summer = new Date(currentYear, 6, 1);
      var winterOffset = winter.getTimezoneOffset();
      var summerOffset = summer.getTimezoneOffset();

      // Local standard timezone offset. Local standard time is not adjusted for
      // daylight savings.  This code uses the fact that getTimezoneOffset
      // returns a greater value during Standard Time versus Daylight Saving
      // Time (DST). Thus it determines the expected output during Standard
      // Time, and it compares whether the output of the given date the same
      // (Standard) or less (DST).
      var stdTimezoneOffset = Math.max(winterOffset, summerOffset);

      // timezone is specified as seconds west of UTC ("The external variable
      // `timezone` shall be set to the difference, in seconds, between
      // Coordinated Universal Time (UTC) and local standard time."), the same
      // as returned by stdTimezoneOffset.
      // See http://pubs.opengroup.org/onlinepubs/009695399/functions/tzset.html
      HEAPU32[((timezone) >> 2)] = stdTimezoneOffset * 60;

      HEAP32[((daylight) >> 2)] = Number(winterOffset != summerOffset);

      var extractZone =
          (timezoneOffset) => {
            // Why inverse sign?
            // Read here
            // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Date/getTimezoneOffset
            var sign = timezoneOffset >= 0 ? "-" : "+";

            var absOffset = Math.abs(timezoneOffset)
            var hours = String(Math.floor(absOffset / 60)).padStart(2, "0");
            var minutes = String(absOffset % 60).padStart(2, "0");

            return `UTC${sign}${hours}${minutes}`;
          }

      var winterName = extractZone(winterOffset);
      var summerName = extractZone(summerOffset);
      assert(winterName);
      assert(summerName);
      assert(lengthBytesUTF8(winterName) <= 16,
             `timezone name truncated to fit in TZNAME_MAX (${winterName})`);
      assert(lengthBytesUTF8(summerName) <= 16,
             `timezone name truncated to fit in TZNAME_MAX (${summerName})`);
      if (summerOffset < winterOffset) {
        // Northern hemisphere
        stringToUTF8(winterName, std_name, 17);
        stringToUTF8(summerName, dst_name, 17);
      } else {
        stringToUTF8(winterName, dst_name, 17);
        stringToUTF8(summerName, std_name, 17);
      }
    };

    var _emscripten_get_now = () => performance.now();

    var _emscripten_date_now = () => Date.now();

    var nowIsMonotonic = 1;

    var checkWasiClock = (clock_id) => clock_id >= 0 && clock_id <= 3;

    function _clock_time_get(clk_id, ignored_precision, ptime) {
      ignored_precision = bigintToI53Checked(ignored_precision);

      if (!checkWasiClock(clk_id)) {
        return 28;
      }
      var now;
      // all wasi clocks but realtime are monotonic
      if (clk_id === 0) {
        now = _emscripten_date_now();
      } else if (nowIsMonotonic) {
        now = _emscripten_get_now();
      } else {
        return 52;
      }
      // "now" is in ms, and wasi times are in ns.
      var nsec = Math.round(now * 1000 * 1000);
      HEAP64[((ptime) >> 3)] = BigInt(nsec);
      return 0;
      ;
    }

    var getHeapMax = () =>
        // Stay one Wasm page short of 4GB: while e.g. Chrome is able to
        // allocate full 4GB Wasm memories, the size will wrap back to 0 bytes
        // in Wasm side for any code that deals with heap sizes, which would
        // require special casing all heap size related code to treat 0
        // specially.
        2147483648;

    var alignMemory = (size, alignment) => {
      assert(alignment, "alignment argument is required");
      return Math.ceil(size / alignment) * alignment;
    };

    var growMemory = (size) => {
      var oldHeapSize = wasmMemory.buffer.byteLength;
      var pages = ((size - oldHeapSize + 65535) / 65536) | 0;
      try {
        // round size grow request up to wasm page size (fixed 64KB per spec)
        wasmMemory.grow(
            pages); // .grow() takes a delta compared to the previous size
        updateMemoryViews();
        return 1 /*success*/;
      } catch (e) {
        err(`growMemory: Attempted to grow heap from ${oldHeapSize} bytes to ${
            size} bytes, but got error: ${e}`);
      }
      // implicit 0 return to save code size (caller will cast "undefined" into
      // 0 anyhow)
    };
    var _emscripten_resize_heap = (requestedSize) => {
      var oldSize = HEAPU8.length;
      // With CAN_ADDRESS_2GB or MEMORY64, pointers are already unsigned.
      requestedSize >>>= 0;
      // With multithreaded builds, races can happen (another thread might
      // increase the size in between), so return a failure, and let the caller
      // retry.
      assert(requestedSize > oldSize);

      // Memory resize rules:
      // 1.  Always increase heap size to at least the requested size, rounded
      // up
      //     to next page multiple.
      // 2a. If MEMORY_GROWTH_LINEAR_STEP == -1, excessively resize the heap
      //     geometrically: increase the heap size according to
      //     MEMORY_GROWTH_GEOMETRIC_STEP factor (default +20%), At most
      //     overreserve by MEMORY_GROWTH_GEOMETRIC_CAP bytes (default 96MB).
      // 2b. If MEMORY_GROWTH_LINEAR_STEP != -1, excessively resize the heap
      //     linearly: increase the heap size by at least
      //     MEMORY_GROWTH_LINEAR_STEP bytes.
      // 3.  Max size for the heap is capped at 2048MB-WASM_PAGE_SIZE, or by
      //     MAXIMUM_MEMORY, or by ASAN limit, depending on which is smallest
      // 4.  If we were unable to allocate as much memory, it may be due to
      //     over-eager decision to excessively reserve due to (3) above.
      //     Hence if an allocation fails, cut down on the amount of excess
      //     growth, in an attempt to succeed to perform a smaller allocation.

      // A limit is set for how much we can grow. We should not exceed that
      // (the wasm binary specifies it, so if we tried, we'd fail anyhow).
      var maxHeapSize = getHeapMax();
      if (requestedSize > maxHeapSize) {
        err(`Cannot enlarge memory, requested ${
            requestedSize} bytes, but the limit is ${maxHeapSize} bytes!`);
        return false;
      }

      // Loop through potential heap size increases. If we attempt a too eager
      // reservation that fails, cut down on the attempted size and reserve a
      // smaller bump instead. (max 3 times, chosen somewhat arbitrarily)
      for (var cutDown = 1; cutDown <= 4; cutDown *= 2) {
        var overGrownHeapSize =
            oldSize * (1 + 0.2 / cutDown); // ensure geometric growth
        // but limit overreserving (default to capping at +96MB overgrowth at
        // most)
        overGrownHeapSize =
            Math.min(overGrownHeapSize, requestedSize + 100663296);

        var newSize = Math.min(
            maxHeapSize,
            alignMemory(Math.max(requestedSize, overGrownHeapSize), 65536));

        var replacement = growMemory(newSize);
        if (replacement) {

          return true;
        }
      }
      err(`Failed to grow the heap from ${oldSize} bytes to ${
          newSize} bytes, not enough memory!`);
      return false;
    };

    var ENV = {};

    var getExecutableName = () => thisProgram || './this.program';
    var getEnvStrings = () => {
      if (!getEnvStrings.strings) {
        // Default values.
        // Browser language detection #8751
        var lang = ((typeof navigator == 'object' && navigator.language) || 'C')
                       .replace('-', '_') +
                   '.UTF-8';
        var env = {
          'USER' : 'web_user',
          'LOGNAME' : 'web_user',
          'PATH' : '/',
          'PWD' : '/',
          'HOME' : '/home/web_user',
          'LANG' : lang,
          '_' : getExecutableName()
        };
        // Apply the user-provided values, if any.
        for (var x in ENV) {
          // x is a key in ENV; if ENV[x] is undefined, that means it was
          // explicitly set to be so. We allow user code to do that to
          // force variables with default values to remain unset.
          if (ENV[x] === undefined)
            delete env[x];
          else
            env[x] = ENV[x];
        }
        var strings = [];
        for (var x in env) {
          strings.push(`${x}=${env[x]}`);
        }
        getEnvStrings.strings = strings;
      }
      return getEnvStrings.strings;
    };

    var _environ_get = (__environ, environ_buf) => {
      var bufSize = 0;
      var envp = 0;
      for (var string of getEnvStrings()) {
        var ptr = environ_buf + bufSize;
        HEAPU32[(((__environ) + (envp)) >> 2)] = ptr;
        bufSize += stringToUTF8(string, ptr, Infinity) + 1;
        envp += 4;
      }
      return 0;
    };

    var _environ_sizes_get = (penviron_count, penviron_buf_size) => {
      var strings = getEnvStrings();
      HEAPU32[((penviron_count) >> 2)] = strings.length;
      var bufSize = 0;
      for (var string of strings) {
        bufSize += lengthBytesUTF8(string) + 1;
      }
      HEAPU32[((penviron_buf_size) >> 2)] = bufSize;
      return 0;
    };

    var keepRuntimeAlive = () => noExitRuntime || runtimeKeepaliveCounter > 0;
    var _proc_exit = (code) => {
      EXITSTATUS = code;
      if (!keepRuntimeAlive()) {
        Module['onExit']?.(code);
        ABORT = true;
      }
      quit_(code, new ExitStatus(code));
    };

    /** @suppress {duplicate } */
    /** @param {boolean|number=} implicit */
    var exitJS = (status, implicit) => {
      EXITSTATUS = status;

      checkUnflushedContent();

      // if exit() was called explicitly, warn the user if the runtime isn't
      // actually being shut down
      if (keepRuntimeAlive() && !implicit) {
        var msg = `program exited (with status: ${
            status}), but keepRuntimeAlive() is set (counter=${
            runtimeKeepaliveCounter}) due to an async operation, so halting execution but not exiting the runtime or preventing further async execution (you can use emscripten_force_exit, if you want to force a true shutdown)`;
        readyPromiseReject?.(msg);
        err(msg);
      }

      _proc_exit(status);
    };
    var _exit = exitJS;

    function _fd_close(fd) {
      try {

        var stream = SYSCALLS.getStreamFromFD(fd);
        FS.close(stream);
        return 0;
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return e.errno;
      }
    }

    function _fd_fdstat_get(fd, pbuf) {
      try {

        var rightsBase = 0;
        var rightsInheriting = 0;
        var flags = 0;
        {
          var stream = SYSCALLS.getStreamFromFD(fd);
          // All character devices are terminals (other things a Linux system
          // would assume is a character device, like the mouse, we have special
          // APIs for).
          var type = stream.tty               ? 2
                     : FS.isDir(stream.mode)  ? 3
                     : FS.isLink(stream.mode) ? 7
                                              : 4;
        }
        HEAP8[pbuf] = type;
        HEAP16[(((pbuf) + (2)) >> 1)] = flags;
        HEAP64[(((pbuf) + (8)) >> 3)] = BigInt(rightsBase);
        HEAP64[(((pbuf) + (16)) >> 3)] = BigInt(rightsInheriting);
        return 0;
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return e.errno;
      }
    }

    /** @param {number=} offset */
    var doReadv = (stream, iov, iovcnt, offset) => {
      var ret = 0;
      for (var i = 0; i < iovcnt; i++) {
        var ptr = HEAPU32[((iov) >> 2)];
        var len = HEAPU32[(((iov) + (4)) >> 2)];
        iov += 8;
        var curr = FS.read(stream, HEAP8, ptr, len, offset);
        if (curr < 0)
          return -1;
        ret += curr;
        if (curr < len)
          break; // nothing more to read
        if (typeof offset != 'undefined') {
          offset += curr;
        }
      }
      return ret;
    };

    function _fd_read(fd, iov, iovcnt, pnum) {
      try {

        var stream = SYSCALLS.getStreamFromFD(fd);
        var num = doReadv(stream, iov, iovcnt);
        HEAPU32[((pnum) >> 2)] = num;
        return 0;
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return e.errno;
      }
    }

    function _fd_seek(fd, offset, whence, newOffset) {
      offset = bigintToI53Checked(offset);

      try {

        if (isNaN(offset))
          return 61;
        var stream = SYSCALLS.getStreamFromFD(fd);
        FS.llseek(stream, offset, whence);
        HEAP64[((newOffset) >> 3)] = BigInt(stream.position);
        if (stream.getdents && offset === 0 && whence === 0)
          stream.getdents = null; // reset readdir state
        return 0;
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return e.errno;
      };
    }

    /** @param {number=} offset */
    var doWritev = (stream, iov, iovcnt, offset) => {
      var ret = 0;
      for (var i = 0; i < iovcnt; i++) {
        var ptr = HEAPU32[((iov) >> 2)];
        var len = HEAPU32[(((iov) + (4)) >> 2)];
        iov += 8;
        var curr = FS.write(stream, HEAP8, ptr, len, offset);
        if (curr < 0)
          return -1;
        ret += curr;
        if (curr < len) {
          // No more space to write.
          break;
        }
        if (typeof offset != 'undefined') {
          offset += curr;
        }
      }
      return ret;
    };

    function _fd_write(fd, iov, iovcnt, pnum) {
      try {

        var stream = SYSCALLS.getStreamFromFD(fd);
        var num = doWritev(stream, iov, iovcnt);
        HEAPU32[((pnum) >> 2)] = num;
        return 0;
      } catch (e) {
        if (typeof FS == 'undefined' || !(e.name === 'ErrnoError'))
          throw e;
        return e.errno;
      }
    }

    var handleException = (e) => {
      // Certain exception types we do not treat as errors since they are used
      // for internal control flow.
      // 1. ExitStatus, which is thrown by exit()
      // 2. "unwind", which is thrown by emscripten_unwind_to_js_event_loop()
      // and others
      //    that wish to return to JS event loop.
      if (e instanceof ExitStatus || e == 'unwind') {
        return EXITSTATUS;
      }
      checkStackCookie();
      if (e instanceof WebAssembly.RuntimeError) {
        if (_emscripten_stack_get_current() <= 0) {
          err('Stack overflow detected.  You can try increasing -sSTACK_SIZE (currently set to 2097152)');
        }
      }
      quit_(1, e);
    };

    var stringToUTF8OnStack = (str) => {
      var size = lengthBytesUTF8(str) + 1;
      var ret = stackAlloc(size);
      stringToUTF8(str, ret, size);
      return ret;
    };

    var AsciiToString = (ptr) => {
      var str = '';
      while (1) {
        var ch = HEAPU8[ptr++];
        if (!ch)
          return str;
        str += String.fromCharCode(ch);
      }
    };

    var FS_createPath = (...args) => FS.createPath(...args);

    var FS_unlink = (...args) => FS.unlink(...args);

    var FS_createLazyFile = (...args) => FS.createLazyFile(...args);

    var FS_createDevice = (...args) => FS.createDevice(...args);

    var incrementExceptionRefcount = (ex) => {
      var ptr = getCppExceptionThrownObjectFromWebAssemblyException(ex);
      ___cxa_increment_exception_refcount(ptr);
    };

    var decrementExceptionRefcount = (ex) => {
      var ptr = getCppExceptionThrownObjectFromWebAssemblyException(ex);
      ___cxa_decrement_exception_refcount(ptr);
    };

    FS.createPreloadedFile = FS_createPreloadedFile;
    FS.preloadFile = FS_preloadFile;
    FS.staticInit();
    ;
    // End JS library code

    // include: postlibrary.js
    // This file is included after the automatically-generated JS library code
    // but before the wasm module is created.

    {

      // Begin ATMODULES hooks
      if (Module['noExitRuntime'])
        noExitRuntime = Module['noExitRuntime'];
      if (Module['preloadPlugins'])
        preloadPlugins = Module['preloadPlugins'];
      if (Module['print'])
        out = Module['print'];
      if (Module['printErr'])
        err = Module['printErr'];
      if (Module['wasmBinary'])
        wasmBinary = Module['wasmBinary'];
      // End ATMODULES hooks

      checkIncomingModuleAPI();

      if (Module['arguments'])
        arguments_ = Module['arguments'];
      if (Module['thisProgram'])
        thisProgram = Module['thisProgram'];

      // Assertions on removed incoming Module JS APIs.
      assert(
          typeof Module['memoryInitializerPrefixURL'] == 'undefined',
          'Module.memoryInitializerPrefixURL option was removed, use Module.locateFile instead');
      assert(
          typeof Module['pthreadMainPrefixURL'] == 'undefined',
          'Module.pthreadMainPrefixURL option was removed, use Module.locateFile instead');
      assert(
          typeof Module['cdInitializerPrefixURL'] == 'undefined',
          'Module.cdInitializerPrefixURL option was removed, use Module.locateFile instead');
      assert(
          typeof Module['filePackagePrefixURL'] == 'undefined',
          'Module.filePackagePrefixURL option was removed, use Module.locateFile instead');
      assert(typeof Module['read'] == 'undefined',
             'Module.read option was removed');
      assert(typeof Module['readAsync'] == 'undefined',
             'Module.readAsync option was removed (modify readAsync in JS)');
      assert(typeof Module['readBinary'] == 'undefined',
             'Module.readBinary option was removed (modify readBinary in JS)');
      assert(
          typeof Module['setWindowTitle'] == 'undefined',
          'Module.setWindowTitle option was removed (modify emscripten_set_window_title in JS)');
      assert(typeof Module['TOTAL_MEMORY'] == 'undefined',
             'Module.TOTAL_MEMORY has been renamed Module.INITIAL_MEMORY');
      assert(
          typeof Module['ENVIRONMENT'] == 'undefined',
          'Module.ENVIRONMENT has been deprecated. To force the environment, use the ENVIRONMENT compile-time option (for example, -sENVIRONMENT=web or -sENVIRONMENT=node)');
      assert(
          typeof Module['STACK_SIZE'] == 'undefined',
          'STACK_SIZE can no longer be set at runtime.  Use -sSTACK_SIZE at link time')
      // If memory is defined in wasm, the user can't provide it, or set
      // INITIAL_MEMORY
      assert(
          typeof Module['wasmMemory'] == 'undefined',
          'Use of `wasmMemory` detected.  Use -sIMPORTED_MEMORY to define wasmMemory externally');
      assert(
          typeof Module['INITIAL_MEMORY'] == 'undefined',
          'Detected runtime INITIAL_MEMORY setting.  Use -sIMPORTED_MEMORY to define wasmMemory dynamically');
    }

    // Begin runtime exports
    Module['addRunDependency'] = addRunDependency;
    Module['removeRunDependency'] = removeRunDependency;
    Module['callMain'] = callMain;
    Module['getValue'] = getValue;
    Module['UTF8ToString'] = UTF8ToString;
    Module['AsciiToString'] = AsciiToString;
    Module['FS_preloadFile'] = FS_preloadFile;
    Module['FS_unlink'] = FS_unlink;
    Module['FS_createPath'] = FS_createPath;
    Module['FS_createDevice'] = FS_createDevice;
    Module['FS'] = FS;
    Module['FS_createDataFile'] = FS_createDataFile;
    Module['FS_createLazyFile'] = FS_createLazyFile;
    Module['WORKERFS'] = WORKERFS;
    Module['PROXYFS'] = PROXYFS;
    var missingLibrarySymbols = [
      'writeI53ToI64',
      'writeI53ToI64Clamped',
      'writeI53ToI64Signaling',
      'writeI53ToU64Clamped',
      'writeI53ToU64Signaling',
      'readI53FromI64',
      'readI53FromU64',
      'convertI32PairToI53',
      'convertI32PairToI53Checked',
      'convertU32PairToI53',
      'getTempRet0',
      'setTempRet0',
      'zeroMemory',
      'withStackSave',
      'inetPton4',
      'inetNtop4',
      'inetPton6',
      'inetNtop6',
      'readSockaddr',
      'writeSockaddr',
      'readEmAsmArgs',
      'jstoi_q',
      'autoResumeAudioContext',
      'getDynCaller',
      'dynCall',
      'runtimeKeepalivePush',
      'runtimeKeepalivePop',
      'callUserCallback',
      'maybeExit',
      'asmjsMangle',
      'HandleAllocator',
      'getNativeTypeSize',
      'addOnInit',
      'addOnPostCtor',
      'addOnPreMain',
      'addOnExit',
      'STACK_SIZE',
      'STACK_ALIGN',
      'POINTER_SIZE',
      'ASSERTIONS',
      'ccall',
      'cwrap',
      'convertJsFunctionToWasm',
      'getEmptyTableSlot',
      'updateTableMap',
      'getFunctionAddress',
      'addFunction',
      'removeFunction',
      'intArrayToString',
      'stringToAscii',
      'UTF16ToString',
      'stringToUTF16',
      'lengthBytesUTF16',
      'UTF32ToString',
      'stringToUTF32',
      'lengthBytesUTF32',
      'stringToNewUTF8',
      'writeArrayToMemory',
      'registerKeyEventCallback',
      'maybeCStringToJsString',
      'findEventTarget',
      'getBoundingClientRect',
      'fillMouseEventData',
      'registerMouseEventCallback',
      'registerWheelEventCallback',
      'registerUiEventCallback',
      'registerFocusEventCallback',
      'fillDeviceOrientationEventData',
      'registerDeviceOrientationEventCallback',
      'fillDeviceMotionEventData',
      'registerDeviceMotionEventCallback',
      'screenOrientation',
      'fillOrientationChangeEventData',
      'registerOrientationChangeEventCallback',
      'fillFullscreenChangeEventData',
      'registerFullscreenChangeEventCallback',
      'JSEvents_requestFullscreen',
      'JSEvents_resizeCanvasForFullscreen',
      'registerRestoreOldStyle',
      'hideEverythingExceptGivenElement',
      'restoreHiddenElements',
      'setLetterbox',
      'softFullscreenResizeWebGLRenderTarget',
      'doRequestFullscreen',
      'fillPointerlockChangeEventData',
      'registerPointerlockChangeEventCallback',
      'registerPointerlockErrorEventCallback',
      'requestPointerLock',
      'fillVisibilityChangeEventData',
      'registerVisibilityChangeEventCallback',
      'registerTouchEventCallback',
      'fillGamepadEventData',
      'registerGamepadEventCallback',
      'registerBeforeUnloadEventCallback',
      'fillBatteryEventData',
      'registerBatteryEventCallback',
      'setCanvasElementSize',
      'getCanvasElementSize',
      'jsStackTrace',
      'getCallstack',
      'convertPCtoSourceLocation',
      'wasiRightsToMuslOFlags',
      'wasiOFlagsToMuslOFlags',
      'safeSetTimeout',
      'setImmediateWrapped',
      'safeRequestAnimationFrame',
      'clearImmediateWrapped',
      'registerPostMainLoop',
      'registerPreMainLoop',
      'getPromise',
      'makePromise',
      'idsToPromises',
      'makePromiseCallback',
      'Browser_asyncPrepareDataCounter',
      'arraySum',
      'addDays',
      'getSocketFromFD',
      'getSocketAddress',
      'FS_mkdirTree',
      '_setNetworkCallback',
      'heapObjectForWebGLType',
      'toTypedArrayIndex',
      'webgl_enable_ANGLE_instanced_arrays',
      'webgl_enable_OES_vertex_array_object',
      'webgl_enable_WEBGL_draw_buffers',
      'webgl_enable_WEBGL_multi_draw',
      'webgl_enable_EXT_polygon_offset_clamp',
      'webgl_enable_EXT_clip_control',
      'webgl_enable_WEBGL_polygon_mode',
      'emscriptenWebGLGet',
      'computeUnpackAlignedImageSize',
      'colorChannelsInGlTextureFormat',
      'emscriptenWebGLGetTexPixelData',
      'emscriptenWebGLGetUniform',
      'webglGetUniformLocation',
      'webglPrepareUniformLocationsBeforeFirstUse',
      'webglGetLeftBracePos',
      'emscriptenWebGLGetVertexAttrib',
      '__glGetActiveAttribOrUniform',
      'writeGLArray',
      'registerWebGlEventCallback',
      'runAndAbortIfError',
      'ALLOC_NORMAL',
      'ALLOC_STACK',
      'allocate',
      'writeStringToMemory',
      'writeAsciiToMemory',
      'demangle',
      'stackTrace',
    ];
    missingLibrarySymbols.forEach(missingLibrarySymbol)

    var unexportedSymbols = [
      'run',
      'out',
      'err',
      'abort',
      'wasmMemory',
      'wasmExports',
      'HEAPF32',
      'HEAPF64',
      'HEAP8',
      'HEAPU8',
      'HEAP16',
      'HEAPU16',
      'HEAP32',
      'HEAPU32',
      'HEAP64',
      'HEAPU64',
      'writeStackCookie',
      'checkStackCookie',
      'INT53_MAX',
      'INT53_MIN',
      'bigintToI53Checked',
      'stackSave',
      'stackRestore',
      'stackAlloc',
      'ptrToString',
      'exitJS',
      'getHeapMax',
      'growMemory',
      'ENV',
      'ERRNO_CODES',
      'strError',
      'DNS',
      'Protocols',
      'Sockets',
      'timers',
      'warnOnce',
      'readEmAsmArgsArray',
      'getExecutableName',
      'handleException',
      'keepRuntimeAlive',
      'asyncLoad',
      'alignMemory',
      'mmapAlloc',
      'wasmTable',
      'getUniqueRunDependency',
      'noExitRuntime',
      'addOnPreRun',
      'addOnPostRun',
      'freeTableIndexes',
      'functionsInTableMap',
      'setValue',
      'PATH',
      'PATH_FS',
      'UTF8Decoder',
      'UTF8ArrayToString',
      'stringToUTF8Array',
      'stringToUTF8',
      'lengthBytesUTF8',
      'intArrayFromString',
      'UTF16Decoder',
      'stringToUTF8OnStack',
      'JSEvents',
      'specialHTMLTargets',
      'findCanvasEventTarget',
      'currentFullscreenStrategy',
      'restoreOldWindowedStyle',
      'UNWIND_CACHE',
      'ExitStatus',
      'getEnvStrings',
      'checkWasiClock',
      'doReadv',
      'doWritev',
      'initRandomFill',
      'randomFill',
      'emSetImmediate',
      'emClearImmediate_deps',
      'emClearImmediate',
      'promiseMap',
      'getExceptionMessageCommon',
      'getCppExceptionTag',
      'getCppExceptionThrownObjectFromWebAssemblyException',
      'Browser',
      'requestFullscreen',
      'requestFullScreen',
      'setCanvasSize',
      'getUserMedia',
      'createContext',
      'getPreloadedImageData__data',
      'wget',
      'MONTH_DAYS_REGULAR',
      'MONTH_DAYS_LEAP',
      'MONTH_DAYS_REGULAR_CUMULATIVE',
      'MONTH_DAYS_LEAP_CUMULATIVE',
      'isLeapYear',
      'ydayFromDate',
      'SYSCALLS',
      'preloadPlugins',
      'FS_createPreloadedFile',
      'FS_modeStringToFlags',
      'FS_getMode',
      'FS_stdin_getChar_buffer',
      'FS_stdin_getChar',
      'FS_readFile',
      'FS_root',
      'FS_mounts',
      'FS_devices',
      'FS_streams',
      'FS_nextInode',
      'FS_nameTable',
      'FS_currentPath',
      'FS_initialized',
      'FS_ignorePermissions',
      'FS_filesystems',
      'FS_syncFSRequests',
      'FS_readFiles',
      'FS_lookupPath',
      'FS_getPath',
      'FS_hashName',
      'FS_hashAddNode',
      'FS_hashRemoveNode',
      'FS_lookupNode',
      'FS_createNode',
      'FS_destroyNode',
      'FS_isRoot',
      'FS_isMountpoint',
      'FS_isFile',
      'FS_isDir',
      'FS_isLink',
      'FS_isChrdev',
      'FS_isBlkdev',
      'FS_isFIFO',
      'FS_isSocket',
      'FS_flagsToPermissionString',
      'FS_nodePermissions',
      'FS_mayLookup',
      'FS_mayCreate',
      'FS_mayDelete',
      'FS_mayOpen',
      'FS_checkOpExists',
      'FS_nextfd',
      'FS_getStreamChecked',
      'FS_getStream',
      'FS_createStream',
      'FS_closeStream',
      'FS_dupStream',
      'FS_doSetAttr',
      'FS_chrdev_stream_ops',
      'FS_major',
      'FS_minor',
      'FS_makedev',
      'FS_registerDevice',
      'FS_getDevice',
      'FS_getMounts',
      'FS_syncfs',
      'FS_mount',
      'FS_unmount',
      'FS_lookup',
      'FS_mknod',
      'FS_statfs',
      'FS_statfsStream',
      'FS_statfsNode',
      'FS_create',
      'FS_mkdir',
      'FS_mkdev',
      'FS_symlink',
      'FS_rename',
      'FS_rmdir',
      'FS_readdir',
      'FS_readlink',
      'FS_stat',
      'FS_fstat',
      'FS_lstat',
      'FS_doChmod',
      'FS_chmod',
      'FS_lchmod',
      'FS_fchmod',
      'FS_doChown',
      'FS_chown',
      'FS_lchown',
      'FS_fchown',
      'FS_doTruncate',
      'FS_truncate',
      'FS_ftruncate',
      'FS_utime',
      'FS_open',
      'FS_close',
      'FS_isClosed',
      'FS_llseek',
      'FS_read',
      'FS_write',
      'FS_mmap',
      'FS_msync',
      'FS_ioctl',
      'FS_writeFile',
      'FS_cwd',
      'FS_chdir',
      'FS_createDefaultDirectories',
      'FS_createDefaultDevices',
      'FS_createSpecialDirectories',
      'FS_createStandardStreams',
      'FS_staticInit',
      'FS_init',
      'FS_quit',
      'FS_findObject',
      'FS_analyzePath',
      'FS_createFile',
      'FS_forceLoadFile',
      'FS_absolutePath',
      'FS_createFolder',
      'FS_createLink',
      'FS_joinPath',
      'FS_mmapAlloc',
      'FS_standardizePath',
      'MEMFS',
      'TTY',
      'PIPEFS',
      'SOCKFS',
      'tempFixedLengthArray',
      'miniTempWebGLFloatBuffers',
      'miniTempWebGLIntBuffers',
      'GL',
      'AL',
      'GLUT',
      'EGL',
      'GLEW',
      'IDBStore',
      'SDL',
      'SDL_gfx',
      'allocateUTF8',
      'allocateUTF8OnStack',
      'print',
      'printErr',
      'jstoi_s',
    ];
    unexportedSymbols.forEach(unexportedRuntimeSymbol);

    // End runtime exports
    // Begin JS library exports
    Module['getExceptionMessage'] = getExceptionMessage;
    Module['incrementExceptionRefcount'] = incrementExceptionRefcount;
    Module['decrementExceptionRefcount'] = decrementExceptionRefcount;
    // End JS library exports

    // end include: postlibrary.js

    function checkIncomingModuleAPI() { ignoredModuleProp('fetchSettings'); }
    function _jsSendStatusUpdateStdout(status_update) {
      postMessage({
        type : "biowasm",
        value : {text : Module.AsciiToString(status_update), type : "print"}
      });
    }
    function _jsSendStatusUpdate(status_update) {
      postMessage({
        type : "biowasm",
        value : {text : Module.AsciiToString(status_update), type : "update"}
      });
    }

    // Imports from the Wasm binary.
    var _free = makeInvalidEarlyAccess('_free');
    var _fflush = makeInvalidEarlyAccess('_fflush');
    var _main = Module['_main'] = makeInvalidEarlyAccess('_main');
    var _emscripten_stack_get_end =
        makeInvalidEarlyAccess('_emscripten_stack_get_end');
    var _emscripten_stack_get_base =
        makeInvalidEarlyAccess('_emscripten_stack_get_base');
    var _strerror = makeInvalidEarlyAccess('_strerror');
    var ___trap = makeInvalidEarlyAccess('___trap');
    var _emscripten_stack_init =
        makeInvalidEarlyAccess('_emscripten_stack_init');
    var _emscripten_stack_get_free =
        makeInvalidEarlyAccess('_emscripten_stack_get_free');
    var __emscripten_stack_restore =
        makeInvalidEarlyAccess('__emscripten_stack_restore');
    var __emscripten_stack_alloc =
        makeInvalidEarlyAccess('__emscripten_stack_alloc');
    var _emscripten_stack_get_current =
        makeInvalidEarlyAccess('_emscripten_stack_get_current');
    var ___cxa_decrement_exception_refcount =
        makeInvalidEarlyAccess('___cxa_decrement_exception_refcount');
    var ___cxa_increment_exception_refcount =
        makeInvalidEarlyAccess('___cxa_increment_exception_refcount');
    var ___thrown_object_from_unwind_exception =
        makeInvalidEarlyAccess('___thrown_object_from_unwind_exception');
    var ___get_exception_message =
        makeInvalidEarlyAccess('___get_exception_message');

    function assignWasmExports(wasmExports) {
      _free = createExportWrapper('free', 1);
      _fflush = createExportWrapper('fflush', 1);
      Module['_main'] = _main = createExportWrapper('__main_argc_argv', 2);
      _emscripten_stack_get_end = wasmExports['emscripten_stack_get_end'];
      _emscripten_stack_get_base = wasmExports['emscripten_stack_get_base'];
      _strerror = createExportWrapper('strerror', 1);
      ___trap = wasmExports['__trap'];
      _emscripten_stack_init = wasmExports['emscripten_stack_init'];
      _emscripten_stack_get_free = wasmExports['emscripten_stack_get_free'];
      __emscripten_stack_restore = wasmExports['_emscripten_stack_restore'];
      __emscripten_stack_alloc = wasmExports['_emscripten_stack_alloc'];
      _emscripten_stack_get_current =
          wasmExports['emscripten_stack_get_current'];
      ___cxa_decrement_exception_refcount =
          createExportWrapper('__cxa_decrement_exception_refcount', 1);
      ___cxa_increment_exception_refcount =
          createExportWrapper('__cxa_increment_exception_refcount', 1);
      ___thrown_object_from_unwind_exception =
          createExportWrapper('__thrown_object_from_unwind_exception', 1);
      ___get_exception_message =
          createExportWrapper('__get_exception_message', 3);
    }
    var ___cpp_exception;
    var wasmImports = {
      /** @export */
      __assert_fail : ___assert_fail,
      /** @export */
      __call_sighandler : ___call_sighandler,
      /** @export */
      __syscall_fcntl64 : ___syscall_fcntl64,
      /** @export */
      __syscall_fstat64 : ___syscall_fstat64,
      /** @export */
      __syscall_getcwd : ___syscall_getcwd,
      /** @export */
      __syscall_ioctl : ___syscall_ioctl,
      /** @export */
      __syscall_lstat64 : ___syscall_lstat64,
      /** @export */
      __syscall_newfstatat : ___syscall_newfstatat,
      /** @export */
      __syscall_openat : ___syscall_openat,
      /** @export */
      __syscall_rmdir : ___syscall_rmdir,
      /** @export */
      __syscall_stat64 : ___syscall_stat64,
      /** @export */
      __syscall_unlinkat : ___syscall_unlinkat,
      /** @export */
      __throw_exception_with_stack_trace : ___throw_exception_with_stack_trace,
      /** @export */
      _abort_js : __abort_js,
      /** @export */
      _emscripten_runtime_keepalive_clear :
          __emscripten_runtime_keepalive_clear,
      /** @export */
      _gmtime_js : __gmtime_js,
      /** @export */
      _jsSendStatusUpdate,
      /** @export */
      _jsSendStatusUpdateStdout,
      /** @export */
      _localtime_js : __localtime_js,
      /** @export */
      _tzset_js : __tzset_js,
      /** @export */
      clock_time_get : _clock_time_get,
      /** @export */
      emscripten_date_now : _emscripten_date_now,
      /** @export */
      emscripten_resize_heap : _emscripten_resize_heap,
      /** @export */
      environ_get : _environ_get,
      /** @export */
      environ_sizes_get : _environ_sizes_get,
      /** @export */
      exit : _exit,
      /** @export */
      fd_close : _fd_close,
      /** @export */
      fd_fdstat_get : _fd_fdstat_get,
      /** @export */
      fd_read : _fd_read,
      /** @export */
      fd_seek : _fd_seek,
      /** @export */
      fd_write : _fd_write,
      /** @export */
      proc_exit : _proc_exit
    };
    var wasmExports = await createWasm();

    // include: postamble.js
    // === Auto-generated postamble setup entry stuff ===

    var calledRun;

    function callMain(args = []) {
      assert(
          runDependencies == 0,
          'cannot call main when async dependencies remain! (listen on Module["onRuntimeInitialized"])');
      assert(typeof onPreRuns === 'undefined' || onPreRuns.length == 0,
             'cannot call main when preRun functions remain to be called');

      var entryFunction = _main;

      args.unshift(thisProgram);

      var argc = args.length;
      var argv = stackAlloc((argc + 1) * 4);
      var argv_ptr = argv;
      args.forEach((arg) => {
        HEAPU32[((argv_ptr) >> 2)] = stringToUTF8OnStack(arg);
        argv_ptr += 4;
      });
      HEAPU32[((argv_ptr) >> 2)] = 0;

      try {

        var ret = entryFunction(argc, argv);

        // if we're not running an evented main loop, it's time to exit
        exitJS(ret, /* implicit = */ true);
        return ret;
      } catch (e) {
        return handleException(e);
      }
    }

    function stackCheckInit() {
      // This is normally called automatically during __wasm_call_ctors but need
      // to get these values before even running any of the ctors so we call it
      // redundantly here.
      _emscripten_stack_init();
      // TODO(sbc): Move writeStackCookie to native to to avoid this.
      writeStackCookie();
    }

    function run(args = arguments_) {

      if (runDependencies > 0) {
        dependenciesFulfilled = run;
        return;
      }

      stackCheckInit();

      preRun();

      // a preRun added a dependency, run will be called later
      if (runDependencies > 0) {
        dependenciesFulfilled = run;
        return;
      }

      function doRun() {
        // run may have just been called through dependencies being fulfilled
        // just in this very frame, or while the async setStatus time below was
        // happening
        assert(!calledRun);
        calledRun = true;
        Module['calledRun'] = true;

        if (ABORT)
          return;

        initRuntime();

        preMain();

        readyPromiseResolve?.(Module);
        Module['onRuntimeInitialized']?.();
        consumedModuleProp('onRuntimeInitialized');

        var noInitialRun = Module['noInitialRun'] || true;
        if (!noInitialRun)
          callMain(args);

        postRun();
      }

      if (Module['setStatus']) {
        Module['setStatus']('Running...');
        setTimeout(() => {
          setTimeout(() => Module['setStatus'](''), 1);
          doRun();
        }, 1);
      } else {
        doRun();
      }
      checkStackCookie();
    }

    function checkUnflushedContent() {
      // Compiler settings do not allow exiting the runtime, so flushing
      // the streams is not possible. but in ASSERTIONS mode we check
      // if there was something to flush, and if so tell the user they
      // should request that the runtime be exitable.
      // Normally we would not even include flush() at all, but in ASSERTIONS
      // builds we do so just for this check, and here we see if there is any
      // content to flush, that is, we check if there would have been
      // something a non-ASSERTIONS build would have not seen.
      // How we flush the streams depends on whether we are in
      // SYSCALLS_REQUIRE_FILESYSTEM=0 mode (which has its own special function
      // for this; otherwise, all the code is inside libc)
      var oldOut = out;
      var oldErr = err;
      var has = false;
      out = err = (x) => { has = true; } try { // it doesn't matter if it fails
        _fflush(0);
        // also flush in the JS FS layer
        ['stdout', 'stderr'].forEach((name) => {
          var info = FS.analyzePath('/dev/' + name);
          if (!info)
            return;
          var stream = info.object;
          var rdev = stream.rdev;
          var tty = TTY.ttys[rdev];
          if (tty?.output?.length) {
            has = true;
          }
        });
      } catch (e) {
      }
      out = oldOut;
      err = oldErr;
      if (has) {
        warnOnce(
            'stdio streams had content in them that was not flushed. you should set EXIT_RUNTIME to 1 (see the Emscripten FAQ), or make sure to emit a newline when you printf etc.');
      }
    }

    function preInit() {
      if (Module['preInit']) {
        if (typeof Module['preInit'] == 'function')
          Module['preInit'] = [ Module['preInit'] ];
        while (Module['preInit'].length > 0) {
          Module['preInit'].shift()();
        }
      }
      consumedModuleProp('preInit');
    }

    preInit();
    run();

    // end include: postamble.js

    // include: postamble_modularize.js
    // In MODULARIZE mode we wrap the generated code in a factory function
    // and return either the Module itself, or a promise of the module.
    //
    // We assign to the `moduleRtn` global here and configure closure to see
    // this as and extern so it won't get minified.

    if (runtimeInitialized) {
      moduleRtn = Module;
    } else {
      // Set up the promise that indicates the Module is initialized
      moduleRtn = new Promise((resolve, reject) => {
        readyPromiseResolve = resolve;
        readyPromiseReject = reject;
      });
    }

    // Assertion for attempting to access module properties on the incoming
    // moduleArg.  In the past we used this object as the prototype of the
    // module and assigned properties to it, but now we return a distinct
    // object.  This keeps the instance private until it is ready (i.e the
    // promise has been resolved).
    for (const prop of Object.keys(Module)) {
      if (!(prop in moduleArg)) {
        Object.defineProperty(moduleArg, prop, {
          configurable : true,
          get() {
            abort(`Access to module property ('${
                prop}') is no longer possible via the module constructor argument; Instead, use the result of the module constructor.`)
          }
        });
      }
    }
    // end include: postamble_modularize.js

    return moduleRtn;
  };
})();

// Export using a UMD style export, or ES6 exports if selected
if (typeof exports === 'object' && typeof module === 'object') {
  module.exports = Module;
  // This default export looks redundant, but it allows TS to import this
  // commonjs style module.
  module.exports.default = Module;
} else if (typeof define === 'function' && define['amd'])
  define([], () => Module);
