


"""
Module name: BioBasis
Introduction: Provide basic tools for my other packages.
Methods: 
    complement : get complement of single base
    reverseComplement : get reverse complement of DNA String
    getSkew : get G-C from the begining to the end
    minimumSkew : get minimum G-C point of the DNA String

"""
module BioBasis

export complement
export reverseComplement
export getSkew
export minimumSkew

"""
Function name: complement
Introduction: get complement of single base
parameters: nucleotide::Char
return comp::Char

Example:
complement('A') -> 'T'
complement('C') -> 'G'
"""
function complement(nucleotide::Char)::Char
    # if comp is 'X', that means some wrong happens
    comp = 'X'
    if nucleotide == 'A'
        comp = 'T'
    elseif nucleotide == 'G'
        comp = 'C'
    elseif nucleotide == 'C'
        comp = 'G'
    elseif nucleotide == 'T'
        comp = 'A'
    end
    return comp
end

"""
Function name: reverseComplement
Introduction: get reverse complement of DNA String
parameters: Text::String
return revCom::String

Example:
complement("CGTA") -> "TACG"
"""
function reverseComplement(Text::String)::String
    revCom = ""
    TextLen = length(Text)
    for i = 1:TextLen
        revCom = complement(Text[i])*revCom
    end
    return revCom
end

"""
Function name: getSkew
Introduction: get G-C from the begining to the end
parameters: Text::String, forcs::Bool = false, if forcs is true, then mark the bgining of the String is 0
return skew::Vector{Int64}

Example:
getSkew("ACCTTGG") -> [0,-1,-2,-2,-2,-1,0]

reference:
Phillip Compeau & Pavel Pevzner. Bioinformatics Algorithms -- An Active Learning Approach, 2nd Edition,  Chapter 1, Page 22.
"""
function getSkew(Text::String, forcs::Bool = false)::Vector{Int64}
    TextLen = length(Text)
    skew = zeros(Int64, TextLen)
    cnt = 0
    for i=1:TextLen
        if Text[i] == 'C'
            cnt -= 1
        elseif Text[i] == 'G'
            cnt += 1
        end
        skew[i] = cnt
    end
    return skew
end

"""
Function name: minimumSkew
Introduction: get minimum G-C point of the DNA String
parameters: Text::String, forcs::Bool = false, if forcs is true, then mark the bgining of the String is 0
return minValuePosition::Vector{Int64}

Example:
getSkew("ACCTTGG") -> [3,4,5]

reference:
Phillip Compeau & Pavel Pevzner. Bioinformatics Algorithms -- An Active Learning Approach, 2nd Edition,  Chapter 1, Page 22.
"""
function minimumSkew(Text::String, forcs::Bool = false)::Vector{Int64}
    TextLen = length(Text)
    skew = getSkew(Text, forcs)
    minValuePosition = Vector{Int64}()
    minValue = minimum(skew)
    for i=1:TextLen
        if skew[i] == minValue
            push!(minValuePosition, i)
        end
    end
    return minValuePosition
end


end 