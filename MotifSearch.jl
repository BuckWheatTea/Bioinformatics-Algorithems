

"""
Module name: MotifSearch
Methods: 

    medianString : Get median string for giving DNAs
    dnasScore : Get DNA Score to see the similarity of the DNAs
    chooseBestAnsOfRandomizedMotifSearch : run RandomizedMotifSearch many times and get optimal motifs
    chooseBestAnsOfGibbsSampler : run simple GibbsSampler many times and get optimal motifs

"""
module MotifSearch

include("BioBasis")
include("BioPatternMatch")

using .BioBasis
using .BioPatternMatch



export medianString
export dnasScore
export chooseBestAnsOfRandomizedMotifSearch
export chooseBestAnsOfGibbsSampler



function motifEnumeration_bruteForce(DnaBox::Vector{String}, patternLength::Int64, maxMisMatch::Int64)::Vector{String} # very low efficient
    allPatterns = getPatterns(patternLength)
    possiblePatterns = Vector{String}()
    for pattern in allPatterns
        jud = true
        for Dna in DnaBox
            if hammingDistanceBetweenPatternAndDna(Dna, pattern) > maxMisMatch
                jud = false
                break
            end
        end
        if jud push!(possiblePatterns, pattern) end
    end
    return possiblePatterns
end

function getProfileMatrix(DnaBox::Vector{String}, laplaceFix::Bool = true)::Array{Float64, 2}
    numOfDna = length(DnaBox)
    lengthOfDna = length(DnaBox[1])
    profileMatrix = zeros(Float64, 4, lengthOfDna)
    for i=1:lengthOfDna
        if !laplaceFix count = Dict{Char, Int64}('A'=> 0, 'C'=>0, 'G'=>0, 'T'=>0)
        elseif laplaceFix count = Dict{Char, Int64}('A'=> 1, 'C'=>1, 'G'=>1, 'T'=>1)
        end
        for j=1:numOfDna
            count[DnaBox[j][i]] += 1
        end
        profileMatrix[1,i] = count['A'] / numOfDna
        profileMatrix[2,i] = count['C'] / numOfDna
        profileMatrix[3,i] = count['G'] / numOfDna
        profileMatrix[4,i] = count['T'] / numOfDna
    end
    return profileMatrix
end

function medianString(DnaBox::Vector{String}, patternLength::Int64)::string
    # Problematic, there should have many answer
    allPatterns = getPatterns(patternLength)
    dis = Inf
    median = Vector{String}()
    for pattern in allPatterns
        thisDis = hammingDistanceBetweenPatternAndDnaBox(DnaBox, pattern)
        if thisDis < dis
            dis = thisDis
            median = pattern
            # push!(median, pattern)
        end
    end
    return median
end

function probablityOfPattern(Pattern::String, ProfileMatrix::Array{Float64, 2})::Float64
    length = size(ProfileMatrix)[2]
    prob = 1
    for i=1:length
        if Pattern[i] == 'A'
            prob*=ProfileMatrix[1, i]
        elseif Pattern[i] == 'C'
            prob*=ProfileMatrix[2,i]
        elseif Pattern[i] == 'G'
            prob*=ProfileMatrix[3,i]
        elseif Pattern[i] == 'T'
            prob*=ProfileMatrix[4,i]
        end
    end
    return prob
end

function mostPropablePattern(ProfileMatrix::Array{Float64, 2}, Text::String)::String
    k = size(ProfileMatrix)[2]
    TextLen = length(Text)
    maxProb = 0
    maxProbPattern = Text[1:k]
    for i=1:TextLen - k + 1
        thisPattern = Text[i:i+k-1]
        thisProb = probablityOfPattern(thisPattern, ProfileMatrix)
        if thisProb > maxProb
            maxProb = thisProb
            maxProbPattern = thisPattern
        end
    end
    return maxProbPattern
end

function dnasScore(DnaBox::Vector{String})::Int64
    numOfDna = length(DnaBox)
    lengthOfDna = length(DnaBox[1])
    score = 0
    for i=1:lengthOfDna
        count = Dict{Char, Int64}('A'=> 0, 'C'=>0, 'G'=>0, 'T'=>0)
        for j=1:numOfDna
            count[DnaBox[j][i]] += 1
        end
        score += numOfDna - findmax(count)[1]
    end
    return score
end



function greedyMotifSearch(DnaBox::Vector{String}, k::Int64, laplaceFix::Bool = true)::Vector{String}
    bestMotifs = [dna[1:k] for dna in DnaBox]
    scoreOfBestMotifs = dnasScore(bestMotifs)
    numOfDna = length(DnaBox)
    lengthOfDna = length(DnaBox[1])
    for i=1:lengthOfDna-k+1
        thisPattern = DnaBox[1][i:i+k-1]
        thisDnaBox = Vector{String}()
        push!(thisDnaBox, thisPattern)
        for j=2:numOfDna
            thisProfileMatrix = getProfileMatrix(thisDnaBox, laplaceFix)
            push!(thisDnaBox, mostPropablePattern(thisProfileMatrix, DnaBox[j]))
        end
        thisScore = dnasScore(thisDnaBox)
        if thisScore < scoreOfBestMotifs
            scoreOfBestMotifs = thisScore
            bestMotifs = thisDnaBox
        end
    end
    return bestMotifs
end

function randomizedMotifSearch(DnaBox::Vector{String}, k::Int64, laplaceFix::Bool = true)
    numOfDna = length(DnaBox)
    lengthOfDna = length(DnaBox[1])
    randomRangeMax = lengthOfDna - k + 1
    bestMotifs = Vector{String}()
    for dna in DnaBox
        randomlyBeginPos = rand(1:randomRangeMax)
        push!(bestMotifs, dna[randomlyBeginPos:randomlyBeginPos+k-1])
    end
    scoreOfBestMotifs = dnasScore(bestMotifs)
    while true
        thisProfileMatrix = getProfileMatrix(bestMotifs, laplaceFix)
        thisMotifs = [mostPropablePattern(thisProfileMatrix, dna) for dna in DnaBox]
        scoreOfThisMotifs = dnasScore(thisMotifs)
        if scoreOfThisMotifs < scoreOfBestMotifs
            bestMotifs = thisMotifs
            scoreOfBestMotifs = scoreOfThisMotifs
        else
            break
        end
    end
    return bestMotifs, scoreOfBestMotifs
end

function chooseBestAnsOfRandomizedMotifSearch(DnaBox::Vector{String}, k::Int64, laplaceFix::Bool = true)
    bestMotifs, bestScore = randomizedMotifSearch(DnaBox, k, laplaceFix)
    for i=1:1000
        thisMotifs, thisScore = randomizedMotifSearch(DnaBox, k, laplaceFix)
        if thisScore < bestScore
            bestScore = thisScore
            bestMotifs = thisMotifs
        end
    end
    return bestMotifs, bestScore
end


function gibbsSampler(DnaBox::Vector{String}, k::Int64, laplaceFix::Bool = true, iterNumber::Int64 = 100)
    numOfDna = length(DnaBox)
    lengthOfDna = length(DnaBox[1])
    randomRangeMax = lengthOfDna - k + 1
    bestMotifs = Vector{String}()
    for dna in DnaBox
        randomlyBeginPos = rand(1:randomRangeMax)
        push!(bestMotifs, dna[randomlyBeginPos:randomlyBeginPos+k-1])
    end
    scoreOfBestMotifs = dnasScore(bestMotifs)
    for j=1:iterNumber
        thisMotifs = copy(bestMotifs)
        indexOfChangeDna = rand(1:numOfDna)
        deleteat!(thisMotifs, indexOfChangeDna)
        thisProfileMatrix = formProfileMatrix(thisMotifs, laplaceFix)
        insert!(thisMotifs, indexOfChangeDna, mostPropablePattern(thisProfileMatrix, DnaBox[indexOfChangeDna]))
        scoreOfThisMotifs = dnasScore(thisMotifs)
        if scoreOfThisMotifs <= scoreOfBestMotifs
            scoreOfBestMotifs = scoreOfThisMotifs
            bestMotifs = thisMotifs
        end
    end
    return bestMotifs, scoreOfBestMotifs
end

function chooseBestAnsOfGibbsSampler(DnaBox::Vector{String}, k::Int64, laplaceFix::Bool = true, iterNumber::Int64 = 100)
    bestMotifs, bestScore = gibbsSampler(DnaBox, k, laplaceFix, iterNumber)
    for i=1:1000
        thisMotifs, thisScore = gibbsSampler(DnaBox, k, laplaceFix, iterNumber)
        if thisScore < bestScore
            bestScore = thisScore
            bestMotifs = thisMotifs
        end
    end
    return bestMotifs, bestScore
end



end