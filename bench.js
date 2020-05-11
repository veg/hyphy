var Benchmark = require('benchmark');
var suite = new Benchmark.Suite;
var execSync = require('child_process').execSync;

// add tests
suite
.add('SmallCodon.bf', function() {
  execSync("/home/runner/work/hyphy/hyphy/HYPHYMP 'tests/hbltests/SimpleOptimizations/SmallCodon.bf'", (error, stdout, stderr) => {
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

// add listeners
.on('cycle', function(event) {
  console.log(String(event.target));
})
.on('complete', function() {
  console.log('Fastest is ' + this.filter('fastest').map('name'));
})
// run async
.run({ 'async': true, 'minSamples': 2});
