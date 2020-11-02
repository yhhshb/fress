## <a name="started"></a>Getting Started
```sh
git clone https://github.com/yhhshb/fress
cd fress && make && make kmc

# Build an SM sketch with approximation factor <epsilon>from a given <kmc database>:
fress sense -i <kmc database> -e <epsilon> -o <output sketch>

# Check errors of a given <sketch>:
fress check -i <kmc database> -d <sketch>
```

## Table of Contents

- [Getting Started](#started)
- [Intro](#uguide)
- [Citing fress](#cite)

## <a name="uguide"></a>Intro

The base mode of operation for fress is to take in input a counting table produced by KMC
and an approximation factor \epsilon in order to generate a SM sketch approximately
associating each k-mer to its frequency.

The main characteristic of the sketches in output is their ability to \epsilon bound the 
total sum of the errors among all k-mers present during construction instead of having
single query guarantees only.

In practice, when the k-mer spectrum of the input data follows a power-law distribution, 
fress can achieve very small errors by relying on the fact that the majority of collisions
will happen between very similar frequencies.

## <a name="cite"></a>Citing fress

If you use minimap2 in your work, please cite:

> Submitted
