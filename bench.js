var Benchmark = require('benchmark');
var _ = require('underscore');
var suite = new Benchmark.Suite;
var execSync = require('child_process').execSync;

function formatNumber(number) {
  number = String(number).split('.');
  return number[0].replace(/(?=(?:\d{3})+$)(?!\b)/g, ',') +
    (number[1] ? '.' + number[1] : '');
}

// add tests
suite
.add('SmallCodon.bf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/SimpleOptimizations/SmallCodon.bf'", (error, stdout, stderr) => {
      if (error) { console.log("Error with job = " + error); return;}
      if (stderr) {return;}
  });
})
.add('IntermediateProtein.bf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/SimpleOptimizations/IntermediateProtein.bf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('FEL.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/FEL.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('SLAC.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/SLAC.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('SLAC-partitioned.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/SLAC-partitioned.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('MEME.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/MEME.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('BUSTED.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/BUSTED.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('BUSTED-SRV.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/BUSTED-SRV.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('RELAX.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/RELAX.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('FUBAR.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/FUBAR.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('BGM.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/BGM.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('CFEL.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/CFEL.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('FADE.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/FADE.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
.add('GARD.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/GARD.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})

.add('ABSREL.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/ABSREL.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})

.add('ABSREL-MH.wbf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/libv3/ABSREL-MH.wbf'", (error, stdout, stderr) => {
    if (error) { console.log("Error with job = " + error); return;}
    if (stderr) {return;}
  });
})
// add listeners
.on('cycle', function(event) {
  let bench = event.target;
  let result = bench.name || (_.isNaN(bench.id) ? bench.id : '<Test #' + bench.id + '>');
  let pm = '\xb1';
  result += ' x ' + formatNumber(bench.hz.toFixed(bench.hz < 100 ? 6 : 0)) + ' ops/sec ' + pm + bench.stats.rme.toFixed(6) + '% (' + bench.stats.sample.length + ' run' + (bench.stats.sample.length == 1 ? '' : 's') + ' sampled)'
  console.log(result);
})
.on('complete', function() {
  console.log('Fastest is ' + this.filter('fastest').map('name'));
})
// run async
.run({ 'async': true });
