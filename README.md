# Synchronization of 3 Christoffel Words
This file contains instructions on how to use the functions to synchronize Christoffel words of three generators.

## Contents:
- orbitize
- christoffelWord
- matrixGenerator
- VI
- transpose
- checkSynchro
- shift
- shiftBy
- allShifts
- shiftOb
- checkFit
- findFits
- merger
- syncWord (Synchronized Word function)
- indexFinder
- synchronizedForm
- seedFinder
- cleanGens
- Synchronize
- synchronizedSeedFunc
- findAllFits
- fitnessSync

### Plotting Functions

- plotCor
- wordPlotter
- InitializePlot2D
- getPoints
- getWord
- plotPoints
- gensPlot2D
- allShiftsGensPlot
- D
- gensReveillesLine

## Details:

- **orbitize**: Finds orbits of any generator.
  - Input: (generator, n = sum of the list of generators)
  - Output: A single list called an orbit of a generator.

- **christoffelWord**: Finds the Christoffel words for a given orbit vector.
  - Input: Orbit for a given generator
  - Output: A single list called the Christoffel word of an orbit.

- **matrixGenerator**: Generates the matrix of orbits and words using the generators.
  - Input: List of generators
  - Output: List of orbits for each generator, list of Christoffel words for each orbit.

- **VI**: Calculates the vertical invariant (d).
  - Input: List of generators
  - Output: Vertical invariant (d).

- **transpose**: Transposes a matrix (list of lists).
  - Input: A matrix (list of lists of equal lengths)
  - Output: A transposed matrix.
    
- **checkSynchro**: check synchronization function: Synchronized iff all sum of cols == 1.
  - Input:  A Christoffel word matrix (list of Christoffel word lists)
  - Output: True or False depending on whether that list is synchronized

- **shift**: Shifts a row by one index.
  - Input: Any list
  - Output: A right shift of the list by one index.
  - Example: input([1,0,0,1,0,0,1]) ----> output([1,1,0,0,1,0,0])

- **shiftBy**: Shifts a list by n indices.
  - Input: Any list
  - Output: A right shift of the list by n indices.

- **allShifts**: generates all possible shifts for a list
  - Input: Any list
  - Output: all possible shifts for that list.

- **shiftOb**: Shifts an orbit by n iterations.
  - Input: An orbit list
  - Output: A right-shifted orbit by n indices.

- **checkFit**: Checks if the sum of columns is either 1 or 0.
  - Input: Christoffel word matrix (or any list of Christoffel words)
  - Output: True or False depending on whether all words fit.

- **findFits**: Finds fits by shifting and checking for fitness at each iteration.
  - Input: Two lists
  - Output: List of all possible fits between them (how many shifts you apply to list 2 for it to fit with list one).

- **merger**: Merges two lists by summing columns.
  - Input: Two lists
  - Output: Lists merged together by summing columns.
  - Example: input([[1,1,1],[0,1,2]]) ----> output([1,2,3])

- **syncWord (Synchronized Word function)**: Finds the synchronized word of a synchronized matrix.
  - Input: Christoffel word matrix
  - Output: Synchronized word, Christoffel Conjugate Matrix.
  - Example: input - [[1,0,0],[0,1,0],[0,0,1]] ----> output - [[[1,0,0],[0,2,0],[0,0,3]], [1,2,3]]

- **indexFinder**: Finds the index of a number in a list. Particularly useful for finding the index of a number in the orbit given the seed.
  - Input: An orbit, a number
  - Output: The index where that number is found.

- **synchronizedForm**: Transforms the orbits using the relative seed (circular permutation starting from a given position).
  - Input: An orbit, an index
  - Output: A shifted orbit such that there is a circular permutation starting from the specified position.
  - Example: synchronizedForm([0,2,4,8,1,3,5,0],3) ---> [8, 1, 3, 5, 0, 2, 4, 8]

- **seedFinder**: Generates a list of seeds.
  - Input: list of generators.
  - output: list of seeds.

- **cleanGens**: Modifies the list of generators (sorted and without zeros).

- **Synchronize**: Uses the seeds to synchronize the Orbit Matrix and words.
  - Input: List of generators, list of seeds
  - Output: Synchronized orbit and Christoffel words matrix.

- **synchronizedSeedFunc**: Checks for synchronization using seeds.
  - Input: A list of generators
  - Output:
    - True/False for synchronization
    - Gens: The generators in order with respect to the output word
    - Words: The synchronized matrix after using the seeds
    - SynchronizedWord: A list consisting of the row numbers where each '1' exists
    - ChristoffelConjugateMatrix: Similar to Words but with '1's replaced by row numbers
    - Orbits: Orbits after applying the seeds
    - Seeds: The seeds (first column of the orbits) used to obtain the synchronized words matrix.

- **findAllFits**: Finds fits in the words matrix in order to synchronize without using seeds.
  - Input: Christoffel word matrix
  - Output: If fits exist, you get a list of all possible fits; if fits don't exist, you get False.

- **fitnessSync**: Checks for synchronization using fits.
  - Input: List of generators
  - Output:
    - True/False: Whether the list of generators can be synchronized or not.
    - CorrectShift: List of shifts applied to achieve synchronization.
    - WorkingMatrix: Synchronized Christoffel words matrix.
    - SynchronizedWord: Synchronized word.
    - ChristoffelConjugateMatrix: Christoffel word matrix with the '1's replaced by row numbers.
    - Generators: List of generators in order of the synchronized Christoffel words matrix.
    - Orbits: Synchronized orbits matrix.
    - Seeds: First column of the synchronized orbits matrix.

### `plotCor`

This function is used to plot a given word based on the starting coordinates.

**Input:**
- `xcor`: Starting x-coordinates
- `ycor`: Starting y-coordinates
- `word`: Christoffel word
- `version`: Version for mapping the words
- `c`: Color of the plot

**Output:**
- Plots the word based on the starting points.

### `wordPlotter`

This function plots a given Christoffel word.

**Input:**
- `word`: Christoffel word
- `c`: Color
- `shift`: A shift to apply to the word
- `version`: Version ('old' or 'new')

**Output:**
- The plot for the word.

### `InitializePlot2D`

This function initializes 2D plots with legends.

**Input:**
- `gens`: List of generators
- `words`: Christoffel word matrix
- `orbits`: Orbital matrix
- `version`: Version ('n' or 'o')
- `c`: Colors used for each row of the Christoffel word matrix
- `shift`: Shift applied

**Output:**
- Plots with legends.

### `getPoints`

This function is used to get all points between the upper and lower bounds.

**Input:**
- `a`, `b`: Slope coordinates
- `mu`: The seed being used

**Output:**
- All the points laying between the bounds.

### `getWord`

This function generates a Christoffel word using the given points.

**Input:**
- `points`: Points laying between bounds

**Output:**
- Christoffel word generated using the points.

### `plotPoints`

This function plots points between bounds.

**Input:**
- `points`: Points laying between bounds
- `a`, `b`: Intercepts of the slope
- `mu`: The seed being used

**Output:**
- A 2D plot for the given points.


### gensPlot2D

**Input:**
- `gens`: List of generators
- `cond`:
  - True: Plotting the synchronized matrix
  - False: Plotting the original unsynchronized Christoffel words matrix
- `version`:
  - Old: Plotting it with 0 meaning moving in the x direction and 1 meaning moving in the y direction
  - New: Similar to old but when we have a 1, map it to 0,1 so we move once in X then once in Y

**Output:**
- A single plot with all the words

### allShiftsGensPlot

**Input:**
- `gens`: List of generators
- `cond`:
  - True: Plotting the synchronized matrix
  - False: Plotting the original unsynchronized Christoffel words matrix
- `version`:
  - Old: Plotting it with 0 meaning moving in the x direction and 1 meaning moving in the y direction
  - New: Similar to old but when we have a 1, map it to 0,1 so we move once in X then once in Y

**Output:**
- n (n = sum of gens) plots with all the words with each plot having shifts applied to it (shifting each plot by i = {1,...,n})

### D

Extracting a word from all the integer points that fall within 2 bounds that can be found using the slope.

**Input:**
- `a`: x coefficient
- `b`: y coefficient
- `mu`: Seed
- `cond`: Whether you want a plot output or not

**Output:**
- The word obtained, plot depending on the condition

### gensReveillesLine

**Input:**
- `gens`: List of generators
- `mu`: Seed applied on the function

**Output:**
- List of words for each generator
- Plots for each word
- Output obtained using the `D` function



## Usage

You can use these functions to analyze and synchronize Christoffel words, explore fitness conditions, and visualize the results.
[Ahmed Abderraouf Menaa]

