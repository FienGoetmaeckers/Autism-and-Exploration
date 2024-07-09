// mean of array
function mean(array) {
  var arraySum = sum(array)
  return arraySum / array.length
}

// sum of array
function sum(array) {
  var total = 0
  for (var i=0; i<array.length; i++) {
      total += array[i]
  }
  return total
}

// repeat array
// repeat complete array a specific number of times
function makeRepeated (arr, repeats) {
  return Array.from({length: repeats}, () => arr).flat()
}