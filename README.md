# Bioinformatics-Algorithms

## Introduction

Using Julia and Python to implement some bioinformatics algorithms mentioned by the book Bioinformatics Algorithms -- An Active Learning Approach, 2nd Edition" by Phillip Compeau &amp; Pavel Pevzner. You can checked detailed introduction in the Julia or Python files. This project are still in progress.

## Important Methods

### BioPatternMatch

- hammingDistance : Computer distance between pattern with same length
- hammingDistanceBetweenPatternAndDna :Computer distance between pattern and DNA String
- hammingDistanceBetweenPatternAndDnaBox :Computer distance between pattern and several DNA Strings.
- approximatePatternMatch : Find similar pattern of given pattern in given DNA String
- approximatePatternCount : To count appear times of given patten in given DNA String, allowed mismatch.
- frequentPatternsWithMismatches : Get most frequent pattern appeared in given DNA String, allowed mismatch

### MotifSearch

- medianString : Get median string for giving DNAs
- dnasScore : Get DNA Score to see the similarity of the DNAs
- chooseBestAnsOfRandomizedMotifSearch : run RandomizedMotifSearch many times and get optimal motifs
- chooseBestAnsOfGibbsSampler : run simple GibbsSampler many times and get optimal motifs

To be continued..

