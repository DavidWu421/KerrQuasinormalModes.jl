using CSV

# Build a lookup table of pre-compiled closures: (s,l,m) => z -> value
# Each formula uses `zm = 1-z` and `zp = 1+z` as locals.
const _SWSH_TABLE = let
    table = Dict{Tuple{Int64,Int64,Int64}, Function}()
    for row in CSV.File(joinpath(@__DIR__, "SpinWeightedSphericalHarmonicsTable.csv"))
        s = Int64(row.s); l = Int64(row.l); m = Int64(row.m)
        # Wrap each formula string into a function body with zm/zp defined
        expr = Meta.parse("""
            let
                local _f = function(z)
                    zm = 1 - z; zp = 1 + z
                    $(row.Func)
                end
                _f
            end
        """)
        table[(s, l, m)] = eval(expr)
    end
    table
end

function sterlings(n)
    return sqrt(2*pi*n)*(n/exp(1))^n
end

function SpinWeightedSphericalCalculation(z, s, l, m)
    if l < max(abs(s), abs(m))
        return Complex(0.0)
    end

    # Fast path: use pre-compiled closure from table
    key = (Int64(s), Int64(l), Int64(m))
    if haskey(_SWSH_TABLE, key)
        return Complex(_SWSH_TABLE[key](z))
    end

    # Fallback: general formula for (s,l,m) outside the table
    if ((l+m > 20) || (l+s > 20) || (l-m > 20) || (l-s > 20))
        term1 = (2*l+1)*(sterlings(l+m)/sterlings(l+s))*(sterlings(l-m)/sterlings(l-s))
    else
        term1 = (2*l+1)*(factorial(l+m)/factorial(l+s))*(factorial(l-m)/factorial(l-s))
    end
    A = sqrt(term1/(4*pi))*(-1)^m

    sinterm = ((1-z)/2)^l
    sumterms = 0.0
    for r = 0:(l-s)
        consts = binomial(l-s,r)*binomial(l+s,r+s-m)*(-1)^(l-r-s)
        cotterm = (sqrt((1+z)/(1-z)))^((2*r+s-m))
        sumterms += consts*cotterm
    end
    return Complex(A*sinterm*sumterms)
end