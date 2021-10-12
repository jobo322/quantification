'use strict';

const optimize = require('./optimizePseudoVoigt');

function optimizePeaks(peakList, x, y, options) {
  let { n = 4, fnType = 'voigt', noiseLevel = 0 } = options;
  let lastIndex = [0];
  let groups = groupPeaks(peakList, n);

  let result = [];
  let sampling, error, opts;

  for (let i = 0; i < groups.length; i++) {
    let peaks = groups[i].group;
    sampling = sampleFunction(
      groups[i].limits[0] - groups[i].limits[1],
      groups[i].limits[0] + groups[i].limits[1],
      x,
      y,
      lastIndex,
    );

    if (sampling[0].length > 5) {
      error = peaks[0].width / 1000;
      opts = {
        LMOptions: [
          3,
          100,
          error,
          error,
          error,
          error * 10,
          error * 10,
          11,
          9,
          1,
        ],
        noiseLevel,
      };
      let optPeaks = [];
      switch (fnType) {
        case 'voigt':
          optPeaks = optimize.optimizePseudoVoigtSum(sampling, peaks, opts);
          break;
        default:
          optPeaks = optimize.optimizePseudoVoigtSum(sampling, peaks, opts);
      }
      let nL = optPeaks[0].length;
      let keys = ['x', 'y', 'width', 'xL'];
      for (let j = 0; j < optPeaks.length; j++) {
        let toPush = {};
        for (let k = 0; k < nL; k++) {
          toPush[keys[k]] = optPeaks[j][k][0];
        }
        toPush.index = peaks[j].index;
        result.push(toPush);
      }
    }
  }

  return result;
}

function groupPeaks(peakList, nL) {
  let group = [];
  let groups = [];
  let i, j;
  let limits = [peakList[0].x, nL * peakList[0].width];
  let upperLimit, lowerLimit;
  // Merge forward
  for (i = 0; i < peakList.length; i++) {
    // If the 2 things overlaps
    if (
      Math.abs(peakList[i].x - limits[0]) <
      nL * peakList[i].width + limits[1]
    ) {
      // Add the peak to the group
      group.push(peakList[i]);
      // Update the group limits
      upperLimit = limits[0] + limits[1];
      if (peakList[i].x + nL * peakList[i].width > upperLimit) {
        upperLimit = peakList[i].x + nL * peakList[i].width;
      }
      lowerLimit = limits[0] - limits[1];
      if (peakList[i].x - nL * peakList[i].width < lowerLimit) {
        lowerLimit = peakList[i].x - nL * peakList[i].width;
      }
      limits = [
        (upperLimit + lowerLimit) / 2,
        Math.abs(upperLimit - lowerLimit) / 2,
      ];
    } else {
      groups.push({ limits: limits, group: group });
      // let optmimalPeak = fitSpectrum(group,limits,spectrum);
      group = [peakList[i]];
      limits = [peakList[i].x, nL * peakList[i].width];
    }
  }
  groups.push({ limits: limits, group: group });
  // Merge backward
  for (i = groups.length - 2; i >= 0; i--) {
    // The groups overlaps
    if (
      Math.abs(groups[i].limits[0] - groups[i + 1].limits[0]) <
      (groups[i].limits[1] + groups[i + 1].limits[1]) / 2
    ) {
      for (j = 0; j < groups[i + 1].group.length; j++) {
        groups[i].group.push(groups[i + 1].group[j]);
      }
      upperLimit = groups[i].limits[0] + groups[i].limits[1];
      if (groups[i + 1].limits[0] + groups[i + 1].limits[1] > upperLimit) {
        upperLimit = groups[i + 1].limits[0] + groups[i + 1].limits[1];
      }
      lowerLimit = groups[i].limits[0] - groups[i].limits[1];
      if (groups[i + 1].limits[0] - groups[i + 1].limits[1] < lowerLimit) {
        lowerLimit = groups[i + 1].limits[0] - groups[i + 1].limits[1];
      }
      // console.log(limits);
      groups[i].limits = [
        (upperLimit + lowerLimit) / 2,
        Math.abs(upperLimit - lowerLimit) / 2,
      ];

      groups.splice(i + 1, 1);
    }
  }
  return groups;
}

function sampleFunction(from, to, x, y, lastIndex) {
  let nbPoints = x.length;
  let sampleX = [];
  let sampleY = [];
  let direction = Math.sign(x[1] - x[0]); // Direction of the derivative
  if (direction === -1) {
    lastIndex[0] = x.length - 1;
  }

  let delta = Math.abs(to - from) / 2;
  let mid = (from + to) / 2;
  let stop = false;
  let index = lastIndex[0];
  while (!stop && index < nbPoints && index >= 0) {
    if (Math.abs(x[index] - mid) <= delta) {
      sampleX.push(x[index]);
      sampleY.push(y[index]);
      index += direction;
    } else {
      // It is outside the range.
      if (Math.sign(mid - x[index]) === 1) {
        // We'll reach the mid going in the current direction
        index += direction;
      } else {
        // There is not more peaks in the current range
        stop = true;
      }
    }
    // console.log(sampleX);
  }
  lastIndex[0] = index;
  return [sampleX, sampleY];
}

module.exports = optimizePeaks;
