## Define the derivative of a qnm function

#Define the derivative of a HeunConfluentRadial function
function ∂r(ψᵣ::HeunConfluentRadial)
    η = ψᵣ.η;
    α = ψᵣ.α;
    ξ = ψᵣ.ξ;
    ζ = ψᵣ.ζ;
    r₊ = ψᵣ.r₊
    r₋ = ψᵣ.r₋
    aₙ = ψᵣ.coeffs
    """Add a sum over different copies of qnm with some
    change"""
    Ψη = HeunConfluentRadial(η-1,α,ξ,ζ,r₊,r₋,aₙ)
    Ψξ = HeunConfluentRadial(η,α,ξ-1,ζ,r₊,r₋,aₙ)
    aₙshift = convert(Vector{eltype(aₙ)},circshift(aₙ,-1))
    aₙshift[end] = zero(eltype(aₙshift))
    nn = 1:length(aₙshift)
    aₙshift = aₙshift .*nn
    aₙshift_static = similar_type(aₙ)(aₙshift)
    Ψaₙ = HeunConfluentRadial(η-1,α+1,ξ,ζ,r₊,r₋,aₙshift_static)
    (im*(η-α))*Ψη + (im*ξ)*Ψξ + ζ*ψᵣ - Ψaₙ
end

#Define the derivation of the qnm function
function ∂r(Ψ::QuasinormalModeFunction)
    s = Ψ.s; l = Ψ.l; m = Ψ.m; n = Ψ.n; a = Ψ.a; ω = Ψ.ω; Alm = Ψ.Alm;
    ψᵣ = Ψ.R
    η = ψᵣ.η;
    α = ψᵣ.α;
    ξ = ψᵣ.ξ;
    ζ = ψᵣ.ζ;
    r₊ = ψᵣ.r₊
    r₋ = ψᵣ.r₋
    aₙ = ψᵣ.coeffs
    """Add a sum over different copies of qnm with some
    change"""
    Ψη = HeunConfluentRadial(η-1,α,ξ,ζ,r₊,r₋,aₙ)
    Ψηf = QuasinormalModeFunction(s,l,m,n,a,ω,Alm,Ψη,Ψ.S)
    Ψξ = HeunConfluentRadial(η,α,ξ-1,ζ,r₊,r₋,aₙ)
    Ψξf = QuasinormalModeFunction(s,l,m,n,a,ω,Alm,Ψξ,Ψ.S)
    aₙshift = convert(Vector{eltype(aₙ)},circshift(aₙ,-1))
    aₙshift[end] = zero(eltype(aₙ))
    nn = 1:length(aₙshift)
    #print("here")
    aₙshift = aₙshift .*nn
    aₙshift_static = similar_type(aₙ)(aₙshift)
    Ψaₙ = HeunConfluentRadial(η-1,α+1,ξ,ζ,r₊,r₋,aₙshift_static)
    Ψaₙf = QuasinormalModeFunction(s,l,m,n,a,ω,Alm,Ψaₙ,Ψ.S)
    (im*(η-α))*Ψηf + (im*ξ)*Ψξf + ζ*Ψ - Ψaₙf
end

function ∂θ(S::SpinWeightedSpheroidal)
    s = S.s; m = S.m; l=S.l; Cllʼ = S.Cllʼ;
    lmin = S.lmin; lmax = S.lmax;
    N = lmax - lmin + 1
    Cllʼp1 = Complex{Float64}[0.0 + 0.0im for _ in 1:N] 
    Cllʼm1 = Complex{Float64}[0.0 + 0.0im for _ in 1:N] 
    for j = 1:N
        lʼ = j+lmin-1;
        term_p1 = (lʼ - s) * (lʼ + s + 1)
        term_m1 = (lʼ + s) * (lʼ - s + 1)

        if term_p1 >= 0
            Cllʼp1[j] = Cllʼ[j] * sqrt(term_p1)
        else
            Cllʼp1[j] = 0.0 + 0.0im
        end

        if term_m1 >= 0
            Cllʼm1[j] = Cllʼ[j] * sqrt(term_m1)
        else
            Cllʼm1[j] = 0.0 + 0.0im
        end
    end
    ψm1 = SpinWeightedSpheroidal(s-1,l,m,Cllʼm1,lmin,lmax)
    ψp1 = SpinWeightedSpheroidal(s+1,l,m,Cllʼp1,lmin, lmax)
    0.5*(ψm1-ψp1)
end

#Define the derivation of the qnm function
function ∂θ(Ψ::QuasinormalModeFunction)
    s = Ψ.s; l = Ψ.l; m = Ψ.m; n = Ψ.n; a = Ψ.a; ω = Ψ.ω; Alm = Ψ.Alm;
    Cllʼ = Ψ.S.Cllʼ; lmin = Ψ.S.lmin; lmax = Ψ.S.lmax;
    N = lmax - lmin + 1
    Cllʼp1 = Complex{Float64}[0.0 + 0.0im for _ in 1:N] 
    Cllʼm1 = Complex{Float64}[0.0 + 0.0im for _ in 1:N] 
    for j = 1:N
        lʼ = j+lmin-1;
        term_p1 = (lʼ - s) * (lʼ + s + 1)
        term_m1 = (lʼ + s) * (lʼ - s + 1)

        if term_p1 >= 0
            Cllʼp1[j] = Cllʼ[j] * sqrt(term_p1)
        else
            Cllʼp1[j] = 0.0 + 0.0im
        end

        if term_m1 >= 0
            Cllʼm1[j] = Cllʼ[j] * sqrt(term_m1)
        else
            Cllʼm1[j] = 0.0 + 0.0im
        end
    end
    
    ψm1 = SpinWeightedSpheroidal(s-1,l,m,Cllʼm1,lmin,lmax)
    ψp1 = SpinWeightedSpheroidal(s+1,l,m,Cllʼp1,lmin, lmax)
    """Add a sum over different copies of qnm with some
    change"""
    Ψm1 = QuasinormalModeFunction(s-1,l,m,n,a,ω,Alm,Ψ.R,ψm1)
    Ψp1 = QuasinormalModeFunction(s+1,l,m,n,a,ω,Alm,Ψ.R,ψp1)
    0.5*(Ψm1-Ψp1)
end

function ∂t(Ψ::QuasinormalModeFunction)
    ω = Ψ.ω;
    -im*ω*Ψ
end

function ∂ϕ(Ψ::QuasinormalModeFunction)
    m = Ψ.m;
    im*m*Ψ
end



function ∂r(Ψ::LinearCombinationOf{T}) where T
    sum(v*∂r(k) for (k,v) in Ψ.dict)
end

function ∂θ(Ψ::LinearCombinationOf{T}) where T
    sum(v ≈ zero(v) ? v : v*∂θ(k) for (k,v) in Ψ.dict)
end

function ∂t(Ψ::LinearCombinationOf{T}) where T
    sum(v ≈ zero(v) ? v : v*∂t(k) for (k,v) in Ψ.dict)
end

function ∂ϕ(Ψ::LinearCombinationOf{T}) where T
    sum(v ≈ zero(v) ? v : v*∂ϕ(k) for (k,v) in Ψ.dict)
end

export ∂r, ∂θ, ∂t, ∂ϕ
