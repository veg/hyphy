var Benchmark = require('benchmark');
var suite = new Benchmark.Suite;
var execSync = require('child_process').execSync;

// add tests
suite
.add('SmallCodon.bf', function() {
  execSync("./HYPHYMP 'tests/hbltests/SimpleOptimizations/SmallCodon.bf'", (error, stdout, stderr) => {
      if (error) { console.log("Error with job = " + error); return;}
      if (stderr) {return;}
  });
})
.add('IntermediateProtein.bf', function() {
  execSync("./HYPHYMP 'tests/hbltests/SimpleOptimizations/IntermediateProtein.bf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('FEL.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/FEL.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('SLAC.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/SLAC.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('SLAC-partitioned.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/SLAC-partitioned.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('MEME.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/MEME.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('BUSTED.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/BUSTED.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('BUSTED-SRV.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/BUSTED-SRV.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('RELAX.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/RELAX.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('FUBAR.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/FUBAR.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('BGM.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/BGM.wbf''", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('CFEL.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/CFEL.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('FADE.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/FADE.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('GARD.wbf', function() {
  execSync("./HYPHYMP 'tests/hbltests/libv3/GARD.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})

// add listeners
.on('cycle', function(event) {
  console.log(String(event.target));
})
.on('complete', function() {
  console.log('Fastest is ' + this.filter('fastest').map('name'));
})
// run async
.run({ 'async': true });
