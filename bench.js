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
  execSync("./HYPHYMP 'tests/hbltests/SimpleOptimizations/SmallCodon.bf'", (error, stdout, stderr) => {
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
.run({ 'async': true, 'minSamples': 2, 'maxSamples':2});
