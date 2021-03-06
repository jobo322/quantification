'use strict';

function getNoiseLevel(y) {
  let mean = 0;

  let stddev = 0;
  let length = y.length;
  for (let i = 0; i < length; ++i) {
    mean += y[i];
  }
  mean /= length;
  let averageDeviations = new Array(length);
  for (let i = 0; i < length; ++i) {
    averageDeviations[i] = Math.abs(y[i] - mean);
  }
  averageDeviations.sort((a, b) => a - b);
  if (length % 2 === 1) {
    stddev = averageDeviations[(length - 1) / 2] / 0.6745;
  } else {
    stddev =
      (0.5 *
        (averageDeviations[length / 2] + averageDeviations[length / 2 - 1])) /
      0.6745;
  }

  return stddev;
}

module.exports = getNoiseLevel;
