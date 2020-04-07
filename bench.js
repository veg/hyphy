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


// add listeners
.on('cycle', function(event) {
  console.log(String(event.target));
})
.on('complete', function() {
  console.log('Fastest is ' + this.filter('fastest').map('name'));
})
// run async
.run({ 'async': true });
