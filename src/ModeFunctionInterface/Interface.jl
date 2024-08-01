using StaticArrays

abstract type Callable end
abstract type CallableAtom <: Callable end
abstract type CallableCombination <: Callable end

### Radial Functions
struct HeunConfluentRadial{T} <: CallableAtom
    η::Complex{Float64}; α::Complex{Float64}; ξ::Complex{Float64}; ζ::Complex{Float64}
    r₊::Float64; r₋::Float64
    coeffs::T
end

function (ψᵣ::HeunConfluentRadial)(r;isconjugate=false,isminus=false)
    η = ψᵣ.η;
    α = ψᵣ.α;
    ξ = ψᵣ.ξ;
    ζ = ψᵣ.ζ;
    r₊ = ψᵣ.r₊;
    r₋ = ψᵣ.r₋
    asymptoticpart = (r₊-r₋)^(α)*(im*(r-r₋))^(η-α)*(im*(r-r₊))^(ξ)*exp(ζ*r)
    x = (r-r₊)/(r-r₋)
    finalsum = Complex(0.0)
    for n in 1:length(ψᵣ.coeffs)
       finalsum += ψᵣ.coeffs[n]*x^(n-1)
    end
    #print(asymptoticpart)
    if isconjugate+isminus==2||isconjugate+isminus==0
        asymptoticpart*finalsum
    elseif isconjugate+isminus==1
        conj(asymptoticpart*finalsum)
    end
end

function ComplexPlot(ψ::HeunConfluentRadial; ztopleft = 0.8 + 1.0im, zbottomright = 2.0 - 0.5im)
    img = portrait(ztopleft, zbottomright, Ψ;
        point_color = cs_d(; colormap=hsv_colors()))
    display(img)
    img
end

### Angular Functions
struct SpinWeightedSpherical
    s::Int64; l::Int64; m::Int64
end

function (Ψ::SpinWeightedSpherical)(z)
    s = Ψ.s; l = Ψ.l; m = Ψ.m;
    if (l >= min(abs(m),abs(s))) & (l >= abs(m))
        return SpinWeightedSphericalCalculation(z,s,l,m)
    end
    return Complex(0.0)
end

function (Ψ::SpinWeightedSpherical)(z,ϕ)
    Ψ(z)*exp(im*Ψ.m*ϕ)
end

struct SpinWeightedSpheroidal{T} <: CallableAtom
    s::Int64; l::Int64; m::Int64
    Cllʼ::T
    lmin::Int64; lmax::Int64
end

function SpinWeightedSpheroidal(s,l,m,Cllʼ)
    lmins = max(abs(s),abs(m));
    lmax = length(Cllʼ) + lmins -1
    SpinWeightedSpheroidal(s,l,m,Cllʼ,lmins, lmax)
end

function SpinWeightedSpheroidalCalculation(z,s,l,m,Cllʼ,lmin,lmax)
    N = lmax - lmin + 1
    val = Complex(0.0)
    for j = 1:N
        lʼ = j+lmin-1;
        val += Cllʼ[j]*SpinWeightedSphericalCalculation(z,s,lʼ,m)
    end
    val
end

function negate_based_on_l(lst, lmin, lmax)
    if abs(lmax-lmin) % 2 == 0
        [i % 2 == 1 ? -x : x for (i, x) in enumerate(lst)]
    else
        [i % 2 == 0 ? -x : x for (i, x) in enumerate(lst)]
    end
end

function (Ψ::SpinWeightedSpheroidal)(z;isconjugate=false,isminus=false)
    s = Ψ.s; l = Ψ.l; m = Ψ.m;
    lmin = Ψ.lmin; lmax = Ψ.lmax; C=Ψ.Cllʼ;
    Cminus=negate_based_on_l(C,lmin,lmax)
    if isconjugate == false
        if isminus==false
            SpinWeightedSpheroidalCalculation(z,s,l,m,Ψ.Cllʼ,lmin,lmax)
        elseif isminus==true
            (-1)^m*conj(SpinWeightedSpheroidalCalculation(z,-s,l,m,Cminus,lmin,lmax))
        end
    elseif isconjugate==true
        if isminus==false
            conj(SpinWeightedSpheroidalCalculation(z,s,l,m,Ψ.Cllʼ,lmin,lmax))
        elseif isminus==true
            (-1)^m*SpinWeightedSpheroidalCalculation(z,-s,l,m,Cminus,lmin,lmax)
        end
    end
end

(Ψ::SpinWeightedSpheroidal)(z,ϕ) = Ψ(z)*exp(im*Ψ.m*ϕ)


function RadialCoefficients(D₀, D₁, D₂, D₃, D₄; N = 250)
    αₙ(n) = (n+1)*(n+D₀)
    βₙ(n) = -2*n^2 + (D₁+2)*n + D₃
    γₙ(n) = (n-1)*(n+D₂-2) + D₄

    an(n) = βₙ(n)/αₙ(n) ; bn(n) = γₙ(n)/αₙ(n);

    ##Set largest rN value
    u1,u2,u3,u4 = rNCoeffs(D₀,D₁,D₂,D₃,D₄)
    rN = 1 + u1*N^(-0.5) + u2*N^(-1) + u3*N^(1.5) + u4*N^(-2)

    ComputeSeriesFromab(an,bn; N=N, rN=rN, PreN=0)
end

function ComputeSeriesFromab(an::Function,bn::Function; N=250, rN = 0.0*im, PreN = 40)
    ##Initialize rn and fn vectors
    rₙ = zeros(Complex{Float64},N+1)
    fₙ = zeros(Complex{Float64},N+1)

    rold = rN
    ##Startup Pass for rₙ
    for n = (N+PreN):-1:(N+1)
        rnew = -bn(n)/(an(n) + rold)
        rold = rnew
    end
    rₙ[N+1] = rold;
    ##Pass for rn
    for n= N:-1:1
        rₙ[n] = -bn(n)/(an(n) + rₙ[n+1])
    end

    fₙ[1] = 1;
    ##Pass for fn
    for n = 2:(N+1)
        fₙ[n] = rₙ[n-1]*fₙ[n-1]
    end
    fₙ
end

### Combining the Mode radial and angular Information
struct QuasinormalModeFunction{T,L} <: CallableAtom
    s::Int64; l::Int64; m::Int64; n::Int64; a::Float64
    ω::Complex{Float64}
    Alm::Complex{Float64}
    R::HeunConfluentRadial{T}
    S::SpinWeightedSpheroidal{L}
end

struct Custom end
function qnmfunction(::typeof(Custom); s=-2,l=2,m=2,n=0,a=0.00, ω = Complex(0.0), Alm = Complex(0.0), Cllʼ = [Complex(0.0)], N=150,isconjugate=false,isminus=false)
    ((ζ,ξ,η),(p,α,γ,δ,σ),(D₀,D₁,D₂,D₃,D₄)) = ParameterTransformations(l,m,s,a,ω,Alm)
    r₊ = 1 + sqrt(1-a^2); r₋ = 1 - sqrt(1-a^2)

    ##Radial WaveFunction
    an = RadialCoefficients(D₀, D₁, D₂, D₃, D₄; N=(N+100))
    an2 = an[1:N];
    aₙ = SVector{length(an2),Complex{Float64}}(an2)
    Ψᵣ = HeunConfluentRadial(η,α,ξ,ζ,r₊,r₋,aₙ)

    ##Angular WaveFunction
    Ψᵪ = SpinWeightedSpheroidal(s,l,m,Cllʼ)

    QuasinormalModeFunction(s,l,m,n,a,ω,Alm,Ψᵣ,Ψᵪ)
end

(Ψ::QuasinormalModeFunction)(r; isconjugate=false,isminus=false) = Ψ.R(r;isconjugate=isconjugate, isminus=isminus) 
(Ψ::QuasinormalModeFunction)(r, z; isconjugate=false, isminus=false) = Ψ.R(r;isconjugate=isconjugate,isminus=isminus) * Ψ.S(z;isconjugate=isconjugate,isminus=isminus)
(Ψ::QuasinormalModeFunction)(r, z, ϕ) =  Ψ.R(r)*Ψ.S(z)*exp(im*Ψ.m*ϕ)
(Ψ::QuasinormalModeFunction)(r, z, ϕ, t) =  Ψ.R(r)*Ψ.S(z)*exp(im*Ψ.m*ϕ)*exp(-im*Ψ.ω*t)

(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r,),Tuple{Number}}) = Ψ.R(x[:r])
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:θ,),Tuple{Number}}) = Ψ.S(cos(x[:θ]))
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:z,),Tuple{Number}}) = Ψ.S(x[:z])
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :θ),Tuple{Number,Number}}) = Ψ.R(x[:r])*Ψ.S(cos(x[:θ]))
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :z),Tuple{Number,Number}}) = Ψ.R(x[:r])*Ψ.S(x[:z])
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :θ, :ϕ),Tuple{Number,Number,Number}}) = Ψ.R(x[:r])*Ψ.S(cos(x[:θ]))*exp(im*Ψ.m*ϕ)
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :z, :ϕ),Tuple{Number,Number,Number}}) = Ψ.R(x[:r])*Ψ.S(x[:z])*exp(im*Ψ.m*ϕ)
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :θ, :ϕ, :t),Tuple{Number,Number,Number,Number}}) = Ψ.R(x[:r])*Ψ.S(cos(x[:θ]))*exp(im*Ψ.m*ϕ)
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :z, :ϕ, :t),Tuple{Number,Number,Number,Number}}) = Ψ.R(x[:r])*Ψ.S(x[:z])*exp(im*Ψ.m*ϕ)

# Construct a Spin sequence that spits out QNM modes
struct ModeSequence
    spin_seq::SpinSequence
end

function (x::ModeSequence)(a::Real; N=150)
    ω, Alm, Cllʼ = x.spin_seq(a)
    @unpack  s,l,m,n = x.spin_seq.mode
    qnmfunction(Custom; s=s,l=l,m=m,n=n,a=a,ω = ω, Alm = Alm, Cllʼ = Cllʼ, N=N)
end

function ModeSequence(;s=-2, l=2, m=2, n=0)
    mode = Mode(s,l,m,n)
    ss = SpinSequence(mode)
    ModeSequence(ss)
end



## Display Functions
import Base.show

function Int2Sub(num,converter)
    if num < 0
       strnum = string(num)[2]
       return "₋"*converter[strnum]
    else
       strnum = string(num)[1]
       return converter[strnum]
    end
end

function Base.show(io::IO, ψ::QuasinormalModeFunction)
    BigDigits = "+-0123456789"
    SmallDigits = "₊₋₀₁₂₃₄₅₆₇₈₉"
    Subdict = Dict(zip(BigDigits,SmallDigits))
    s = Int2Sub(ψ.s,Subdict)
    l = Int2Sub(ψ.l,Subdict)
    m = Int2Sub(ψ.m,Subdict)
    n = Int2Sub(ψ.n,Subdict)
    print(io,s*"Ψ"*l*m*n)
end
