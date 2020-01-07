


"""
Module name: PatternMatch
Methods:

    patternCount : To count appear times of given patten in given DNA String, mismatch is not allowed
    patternInGenome : Get all start position of given pattern in given DNA String
    frequentWords : Get most frequent pattern appeared in given DNA String
    hammingDistance : Computer distance between pattern with same length
    approximatePatternMatch : Find similar pattern of given pattern in given DNA String
    approximatePatternCount : To count appear times of given patten in given DNA String, allowed mismatch.
    getPatterns : Enumerate all patterns with given length
    frequentPatternsWithMismatches : Get most frequent pattern appeared in given DNA String, allowed mismatch
    frequentPatternsWithMismatchesAndReverseComplements :
    hammingDistanceBetweenPatternAndDna :Computer distance between pattern and DNA String
    hammingDistanceBetweenPatternAndDnaBox :Computer distance between pattern and several DNA Strings.

"""
module BioPatternMatch

include("BioBasis.jl")
using .BioBasis


export patternCount
export patternInGenome
export frequentWords
export hammingDistance
export approximatePatternMatch
export approximatePatternCount
export getPatterns
export frequentPatternsWithMismatches
export frequentPatternsWithMismatchesAndReverseComplements
export hammingDistanceBetweenPatternAndDna
export hammingDistanceBetweenPatternAndDnaBox


"""
Function name : patternCount

Introduction : To count appear times of given patten in given DNA String, mismatch is not allowed

parameters: nucleotide::Text::String, Pattern::String

return : count::Int64

Example:

"""
function patternCount(Text::String, Pattern::String)::Int64
    count = 0
    TextLen = length(Text)
    PatternLen = length(Pattern)
    for i = 1:TextLen-PatternLen+1
        if Text[i:i+PatternLen-1] == Pattern
            count+=1
        end
    end
    return count
end

"""
Function name: patternInGenome
Introduction: Get all start position of given pattern in given DNA String
parameters: Text::String, pattern::String, forcs::Bool=false
return pos::Vector{Int64}

Example:

"""
function patternInGenome(Text::String, pattern::String, forcs::Bool=false)::Vector{Int64}
    # forcs: show the vector for Computer scientist, the index will begin with 0
    pos = Vector{Int64}()
    TextLen = length(Text)
    patternLen = length(pattern)
    for i = 1:TextLen-patternLen+1
        if Text[i:i+patternLen-1] == pattern
            push!(pos, i)
        end
    end
    if forcs
        pos = [i-1 for i in pos]
    end
    return pos
end


"""
Function name: frequentWords
Introduction: Get most frequent pattern appeared in given DNA String
parameters: Text::String, k::Int64
return mostFrequentPatterns::Vector{String}

Example:

"""
function frequentWords(Text::String, k::Int64)::Vector{String}
    frequentPatterns = Dict{String, Int64}()
    TextLen = length(Text)
    max = 0
    for i = 1:TextLen-k+1
        key = Text[i:i+k-1]
        if haskey(frequentPatterns, key)
            frequentPatterns[key] += 1
        else
            frequentPatterns[key] = 1
        end
        if frequentPatterns[key] > max
            max = frequentPatterns[key]
        end
    end
    mostFrequentPatterns = Vector{String}()
    for i in frequentPatterns
        if i.second == max
            push!(mostFrequentPatterns, i.first)
        end
    end
    return mostFrequentPatterns
end

# function complement(nucleotide::Char)::Char
#     # if comp is 'X', that means some wrong happens
#     comp = 'X'
#     if nucleotide == 'A'
#         comp = 'T'
#     elseif nucleotide == 'G'
#         comp = 'C'
#     elseif nucleotide == 'C'
#         comp = 'G'
#     elseif nucleotide == 'T'
#         comp = 'A'
#     end
#     return comp
# end

# function reverseComplement(Text::String)::String
#     revCom = ""
#     TextLen = length(Text)
#     for i = 1:TextLen
#         revCom = complement(Text[i])*revCom
#     end
#     return revCom
# end



# function isClumpApperaTtimes(posVec::Vector{Int64}, L, t)::Bool
#     M = length(posVec)
#     if M < t return false end
#     if posVec[M] - posVec[1] <= L return true
#     else return isClumpApperaTtimes(posVec[1:M-1], L, t) | isClumpApperaTtimes(posVec[2:M], L, t)
#     end
# end

# function getClumps(Text::String, k, L, t)::Vector{String}
#     TextLen = length(Text)
#     fragCountDict = Dict{String, Vector{Int64}}()
#     clumps = Vector{String}()
#     for i=1:TextLen-k+1
#         thisFrag = Text[i:i+k-1]
#         if haskey(fragCountDict, thisFrag) # may low efficient
#             push!(fragCountDict[thisFrag], i)
#         else
#             fragCountDict[thisFrag] = [i,]
#         end
#     end
#     for frag in keys(fragCountDict)
#         if isClumpApperaTtimes(fragCountDict[frag], L-k, t)
#             push!(clumps, frag)
#         end
#     end
#     return clumps
# end

# function getSkew(Text::String, forcs::Bool = false)::Vector{Int64}
#     TextLen = length(Text)
#     skew = zeros(Int64, TextLen)
#     cnt = 0
#     for i=1:TextLen
#         if Text[i] == 'C'
#             cnt -= 1
#         elseif Text[i] == 'G'
#             cnt += 1
#         end
#         skew[i] = cnt
#     end
#     return skew
# end

# function minimumSkew(Text::String, forcs::Bool = false)::Vector{Int64}
#     TextLen = length(Text)
#     skew = getSkew(Text, forcs)
#     minValuePosition = Vector{Int64}()
#     minValue = minimum(skew)
#     for i=1:TextLen
#         if skew[i] == minValue
#             push!(minValuePosition, i)
#         end
#     end
#     return minValuePosition
# end

"""
Function name: hammingDistance
Introduction: Get most frequent pattern appeared in given DNA String
parameters: Text1::String, Text2::String
return hamDis::Int64

Example:

"""
function hammingDistance(Text1::String, Text2::String)::Int64
    TextLen1 = length(Text1)
    TextLen2 = length(Text2)
    if TextLen1 != TextLen2
        println("Wrong! The two text do not have the same length!")
        return -1
    end
    hamDis = 0
    for i=1:TextLen1
        if Text1[i] != Text2[i]
            hamDis+=1
        end
    end
    return hamDis
end


function approximatePatternMatch(Text::String, pattern::String, d::Int64, forcs::Bool = false)::Vector{Int64}
    TextLen = length(Text)
    patternLen = length(pattern)
    posVec = Vector{Int64}()
    for i=1:TextLen - patternLen + 1
        if hammingDistance(Text[i:i+patternLen-1], pattern) <= d
            if forcs push!(posVec, i-1)
            else push!(posVec, i)
            end
        end
    end
    return posVec
end

function approximatePatternCount(Text::String, pattern::String, d::Int64)::Int64
    TextLen = length(Text)
    patternLen = length(pattern)
    cnt = 0
    for i=1:TextLen - patternLen + 1
        if hammingDistance(Text[i:i+patternLen-1], pattern) <= d
            cnt+=1
        end
    end
    return cnt
end


function getPatterns(n::Int64, lastVector::Vector{String} = Vector{String}())::Vector{String} #low efficient
    if n == 0
        return lastVector
    end
    if length(lastVector) == 0
        thisVec = Vector{String}()
        push!(thisVec, "A")
        push!(thisVec, "T")
        push!(thisVec, "C")
        push!(thisVec, "G")
        return getPatterns(n-1, thisVec)
    else
        temVec = ["A", "T", "C", "G"]
        thisVec = Vector{String}()
        for i in temVec
            for j in lastVector
                push!(thisVec, i*j)
            end
        end
        return getPatterns(n-1, thisVec)
    end
end

function frequentPatternsWithMismatches(Text::String, n::Int64, d::Int64)::Vector{String} # very low efficient
    patterns = getPatterns(n)
    patternsCount = Dict{String, Int64}()
    ansPattern = Vector{String}()
    maxValue = 0
    for pattern in patterns
        cnt = approximatePatternCount(Text, pattern, d)
        patternsCount[pattern] = cnt
        if maxValue < cnt
            maxValue = cnt
        end
    end
    for pair in patternsCount
        if pair.second == maxValue
            push!(ansPattern, pair.first)
        end
    end
    return ansPattern
end

function frequentPatternsWithMismatchesAndReverseComplements(Text::String, n::Int64, d::Int64)::Vector{String} # very low efficient
    patterns = getPatterns(n)
    patternsCount = Dict{String, Int64}()
    ansPattern = Vector{String}()
    maxValue = 0
    for pattern in patterns
        revCom = reverseComplement(pattern)
        cnt1 = approximatePatternCount(Text, pattern, d)
        cnt2 = approximatePatternCount(Text, revCom, d)
        cnt = cnt1 + cnt2
        patternsCount[pattern] = cnt
        if maxValue < cnt
            maxValue = cnt
        end
    end
    for pair in patternsCount
        if pair.second == maxValue
            push!(ansPattern, pair.first)
        end
    end
    return ansPattern
end

function hammingDistanceBetweenPatternAndDna(Text::String, pattern::String)::Int64
    TextLen = length(Text)
    patternLen = length(pattern)
    dis = Inf
    for i=1:TextLen - patternLen + 1
        thisDis = hammingDistance(Text[i:i+patternLen-1], pattern)
        if thisDis < dis  dis = thisDis end
    end
    return dis
end

function hammingDistanceBetweenPatternAndDnaBox(DnaBox::Vector{String}, pattern::String)::Int64
    totalDis = 0
    for Dna in DnaBox
        totalDis+=hammingDistanceBetweenPatternAndDna(Dna, pattern)
    end
    return totalDis
end






end

