var Benchmark = require('benchmark');
const { exec } = require("child_process");
var suite = new Benchmark.Suite;


// add tests
suite.add('SmallCodon.bf', function() {
  exec("./HYPHYMP 'tests/hbltests/SimpleOptimizations/SmallCodon.bf'", (error, stdout, stderr) => {
      if (error) {return;}
      if (stderr) {return;}
  });

})
.add('IntermediateProtein.bf', function() {
  exec("./HYPHYMP 'tests/hbltests/SimpleOptimizations/IntermediateProtein.bf'", (error, stdout, stderr) => {
    if (error) {return;}
    if (stderr) {return;}
  });
})


.add('FEL.wbf', function() {
  exec("./HYPHYMP 'tests/hbltests/libv3/FEL.wbf'", (error, stdout, stderr) => {
    if (error) {return;}
    if (stderr) {return;}
  });
})


.add('SLAC.wbf', function() {
  exec("./HYPHYMP 'tests/hbltests/libv3/SLAC.wbf'", (error, stdout, stderr) => {
    if (error) {return;}
    if (stderr) {return;}
  });
})


.add('SLAC-partitioned.wbf', function() {
  exec("./HYPHYMP 'tests/hbltests/libv3/SLAC-partitioned.wbf'", (error, stdout, stderr) => {
    if (error) {return;}
    if (stderr) {return;}
  });
})


.add('MEME.wbf', function() {
  exec("./HYPHYMP 'tests/hbltests/libv3/MEME.wbf'", (error, stdout, stderr) => {
    if (error) {return;}
    if (stderr) {return;}
  });
})

.add('BUSTED.wbf', function() {
  exec("./HYPHYMP 'tests/hbltests/libv3/BUSTED.wbf'", (error, stdout, stderr) => {
    if (error) {return;}
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
