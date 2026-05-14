#Functions for the \scriptD_0 and \scriptL^\dagger_0 derivatives from Chandra.

function 𝒟0(Ψ::QuasinormalModeFunction)
    s = Ψ.s; l = Ψ.l; m = Ψ.m; n = Ψ.n; a = Ψ.a; ω = Ψ.ω; Alm = Ψ.Alm;
    ψᵣ = Ψ.R
    η = ψᵣ.η;
    α = ψᵣ.α;
    ξ = ψᵣ.ξ;
    ζ = ψᵣ.ζ;
    r₊ = ψᵣ.r₊
    r₋ = ψᵣ.r₋
    aₙ = ψᵣ.coeffs
    is_conjugate=Ψ.is_conjugate
    is_minus=Ψ.is_minus
    """Add a sum over different copies of qnm with some
    change"""
    # 𝒟0[Ψ]=∂r(Ψ)-iωa(r-r₊)^(-1)(r-r₋)^(-1)+iam(r-r₊)^(-1)(r-r₋)^(-1)
    # -iω(r-r₊)^(-1)(r-r₋)z-iωr₊(r-r₊)^(-1)(r-r₋)^(-1)-iω(r-r₋)(r₊+r₋)\)
    ∂rΨf = ∂r(Ψ)

    Ψ1= HeunConfluentRadial(η-1,α+1,ξ,ζ,r₊,r₋,aₙ,is_conjugate,is_minus)
    Ψ1f = QuasinormalModeFunction(s,l,m,n,a,ω,Alm,Ψ1,Ψ.S,is_conjugate,is_minus)

    aₙ2 = vcat(zero(eltype(aₙ)), aₙ)
    Ψ2= HeunConfluentRadial(η-1,α-1,ξ,ζ,r₊,r₋,aₙ2,is_conjugate,is_minus)
    Ψ2f = QuasinormalModeFunction(s,l,m,n,a,ω,Alm,Ψ2,Ψ.S,is_conjugate,is_minus)

    Ψ3= HeunConfluentRadial(η,α+1,ξ,ζ,r₊,r₋,aₙ,is_conjugate,is_minus)
    Ψ3f = QuasinormalModeFunction(s,l,m,n,a,ω,Alm,Ψ3,Ψ.S,is_conjugate,is_minus)

    if (is_minus==false && is_conjugate==false) || (is_minus==true && is_conjugate==true)
        # ∂rΨ + (im*ω*a^2/(r₊-r₋))*Ψ1 - (im*a*m/(r₊-r₋))*Ψ2 - im*ω*(r₊-r₋)*Ψ3 + (im*ω*r₊^2/(r₊-r₋))*Ψ4 - (im*ω*(r₊+r₋)/(r₊-r₋))*Ψ5
        ∂rΨf + im*(ω*a^2 - a*m + ω*r₊^2)*Ψ1f - im*ω*(r₊-r₋)*Ψ2f - (im*ω*(r₊+r₋)/(r₊-r₋))*Ψ3f
    elseif (is_minus==true && is_conjugate==false) || (is_minus==false && is_conjugate==true)
        # ∂rΨ + (im*ω*a^2/(r₊-r₋))*Ψ1 - (im*a*m/(r₊-r₋))*Ψ2 - im*ω*(r₊-r₋)*Ψ3 + (im*ω*r₊^2/(r₊-r₋))*Ψ4 + (im*ω*(r₊+r₋)/(r₊-r₋))*Ψ5
        ∂rΨf + im*(ω*a^2 - a*m + ω*r₊^2)*Ψ1f - im*ω*(r₊-r₋)*Ψ2f + (im*ω*(r₊+r₋)/(r₊-r₋))*Ψ3f
    end

end

# function ℒ†0

# end

# export 𝒟0, ℒ†0
export 𝒟0